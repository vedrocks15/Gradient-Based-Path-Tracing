#pragma once

// The simplest volumetric renderer:
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implement this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    std::optional<PathVertex> isect_ = intersect(scene, ray, ray_diff);
    if (isect_) {
        PathVertex isect = *isect_;
        int medium_id = isect.exterior_medium_id;
        const Medium &medium = scene.media[medium_id];
        Spectrum sigma_a = get_sigma_a(medium, isect.position);
        Real t = distance(ray.org, isect.position);
        Spectrum transmittance = exp(-sigma_a * t);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[isect.shape_id])) {
            Le = emission(isect, -ray.dir, scene);
        }
        return transmittance * Le;
    }
    return make_zero_spectrum();
}

// The second simplest volumetric renderer:
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implement this!
    
    
    
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    std::optional<PathVertex> isect_ = intersect(scene, ray, ray_diff);

    const Medium &medium = scene.media[scene.camera.medium_id];
    Real sigma_s = get_sigma_s(medium, ray.org)[0];  // Homogeneous and monochromatic medium
    Real sigma_a = get_sigma_a(medium, ray.org)[0];  // Homogeneous and monochromatic medium
    Real sigma_t = sigma_s + sigma_a;                // Monochromatic medium

    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1 - u) / sigma_t;

    if (!isect_ || t < distance(ray.org, isect_->position)) {
        Real trans_pdf = exp(-sigma_t * t) * sigma_t;
        Real transmittance = exp(-sigma_t * t);
        Vector3 p = ray.org + t * ray.dir;

        // Compute equation (7)
        // Sample a point on the light source
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
        Vector3 dir_light = normalize(point_on_light.position - p);
        Real L_s1_pdf = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, p, scene);

        Spectrum L_s1_estimate;

        // Test if occluded
        Ray shadow_ray{p, dir_light,
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                           distance(point_on_light.position, p)};
        
        if (occluded(scene, shadow_ray)) {
            L_s1_estimate = make_zero_spectrum();
        } else {
            // Evaluate emission
            Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);

            // Evaluate phase function
            PhaseFunction phase = get_phase_function(medium);
            Spectrum phase_val = eval(phase, -ray.dir, dir_light);

            Real dist = distance(p, point_on_light.position);
            Real cos_light = max(-dot(dir_light, point_on_light.normal), Real(0));
            L_s1_estimate = phase_val * Le * exp(-sigma_t * dist) * cos_light / (dist * dist);
        }

        return (transmittance / trans_pdf) * sigma_s * (L_s1_estimate / L_s1_pdf);
    } else {
        Real trans_pdf = exp(-sigma_t * distance(ray.org, isect_->position));
        Real transmittance = exp(-sigma_t * distance(ray.org, isect_->position));
        Spectrum Le = make_zero_spectrum();
        if (isect_ && is_light(scene.shapes[isect_->shape_id])) {
            Le = emission(*isect_, -ray.dir, scene);
        }

        return (transmittance / trans_pdf) * Le;
    }
}

int update_medium_id(const PathVertex &isect, Ray &ray, int current_medium_id) {
    int new_medium_id = current_medium_id;
    if (isect.interior_medium_id != isect.exterior_medium_id) {
        // At medium transition, update medium
        if (dot(ray.dir, isect.geometric_normal) > 0) {
            new_medium_id = isect.exterior_medium_id;
        } else {
            new_medium_id = isect.interior_medium_id;
        }
    }
    return new_medium_id;
}

// The third volumetric renderer (not so simple anymore):
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implement this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = make_const_spectrum(1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    while (true) {
        bool scatter = false;
        std::optional<PathVertex> isect_ = intersect(scene, ray, ray_diff);
        Real transmittance = 1;  // Monochromatic medium
        Real trans_pdf = 1;
        if (current_medium_id != -1) {  // Not in vacuum
            const Medium &medium = scene.media[current_medium_id];
            Real sigma_s = get_sigma_s(medium, ray.org)[0];  // Homogeneous and monochromatic medium
            Real sigma_a = get_sigma_a(medium, ray.org)[0];  // Homogeneous and monochromatic medium
            Real sigma_t = sigma_s + sigma_a;                // Monochromatic medium

            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;

            if (!isect_ || t < distance(ray.org, isect_->position)) {
                scatter = true;
                transmittance = exp(-sigma_t * t);
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                ray.org += t * ray.dir;
            } else {
                transmittance = exp(-sigma_t * distance(ray.org, isect_->position));
                trans_pdf = exp(-sigma_t * distance(ray.org, isect_->position));
                ray.org = isect_->position + ray.dir * get_intersection_epsilon(scene);
            }
        }

        current_path_throughput *= (transmittance / trans_pdf);

        if (!scatter && isect_ && is_light(scene.shapes[isect_->shape_id])) {
            // Reach an emissive surface
            radiance += current_path_throughput * emission(*isect_, -ray.dir, scene);
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) {
            // Reach maximum bounces
            break;
        }

        if (!scatter && isect_ && isect_->material_id == -1) {
            // Index-matching surface, skip through it
            if (isect_->interior_medium_id != isect_->exterior_medium_id) {
                current_medium_id = update_medium_id(*isect_, ray, current_medium_id);
            }
            bounces++;
            continue;
        }

        // Sample next direction and update path throughput
        if (scatter) {
            const Medium &medium = scene.media[current_medium_id];
            Real sigma_s = get_sigma_s(medium, ray.org)[0];  // Homogeneous and monochromatic medium
            PhaseFunction phase = get_phase_function(medium);
            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase, -ray.dir, rnd_param);
            current_path_throughput *= eval(phase, -ray.dir, *next_dir_) /
                                       pdf_sample_phase(phase, -ray.dir, *next_dir_) * sigma_s;
            ray.dir = *next_dir_;
        } else {
            // Hit a surface -- don't need to deal with this yet
            break;
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }

    return radiance;
}

inline Spectrum next_event_estimation_scatter(const Scene &scene, Vector3 p, Vector3 dir_view,
                                              int current_medium_id, int bounces, pcg32_state &rng) {
    // Sample a point on the light source
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 p_prime = point_on_light.position;
    Vector3 dir_light = normalize(p_prime - p);
    Real dist_light = distance(p, p_prime);
    Vector3 init_p = p;

    
    Spectrum T_light = make_const_spectrum(1);  // Transmittance to light
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Real p_trans_dir = 1;
    while (true) {
        Ray shadow_ray{p, dir_light,
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                           distance(p_prime, p)};
        std::optional<PathVertex> isect_ = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime);
        if (isect_) {
            next_t = distance(p, isect_->position);
        }
        if (shadow_medium_id != -1) {
            const Medium &medium = scene.media[shadow_medium_id];
            Real sigma_s = get_sigma_s(medium, p)[0];  // Homogeneous and monochromatic medium
            Real sigma_a = get_sigma_a(medium, p)[0];  // Homogeneous and monochromatic medium
            Real sigma_t = sigma_s + sigma_a;          // Monochromatic medium
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }

        if (!isect_) {
            break;
        } else {
            if (isect_->material_id >= 0) {
                return make_zero_spectrum();
            }
            shadow_bounces++;
            if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                return make_zero_spectrum();
            }
        }

        shadow_medium_id = update_medium_id(*isect_, shadow_ray, shadow_medium_id);
        p = p + next_t * shadow_ray.dir;
    }

    if (max(T_light) > 0) {
        Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) / (dist_light * dist_light);
        PhaseFunction phase = get_phase_function(scene.media[current_medium_id]);
        Spectrum f = eval(phase, dir_view, dir_light);
        Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
        Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, init_p, scene);
        Spectrum contrib = T_light * G * f * L / pdf_nee;

        Real pdf_phase = pdf_sample_phase(phase, dir_view, dir_light) * G * p_trans_dir;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        return w * contrib;
    }

    return make_zero_spectrum();
}

// The fourth volumetric renderer:
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implement this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = make_const_spectrum(1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    Real dir_pdf = 0;
    Vector3 nee_p_cache{0, 0, 0};
    Real multi_trans_pdf = 1;
    bool never_scatter = true;
    while (true) {
        bool scatter = false;
        std::optional<PathVertex> isect_ = intersect(scene, ray, ray_diff);
        Real transmittance = 1;  // Monochromatic medium
        Real trans_pdf = 1;
        if (current_medium_id != -1) {  // Not in vacuum
            const Medium &medium = scene.media[current_medium_id];
            Real sigma_s = get_sigma_s(medium, ray.org)[0];  // Homogeneous and monochromatic medium
            Real sigma_a = get_sigma_a(medium, ray.org)[0];  // Homogeneous and monochromatic medium
            Real sigma_t = sigma_s + sigma_a;                // Monochromatic medium

            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;

            if (!isect_ || t < distance(ray.org, isect_->position)) {
                scatter = true;
                transmittance = exp(-sigma_t * t);
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                ray.org += t * ray.dir;
            } else {
                transmittance = exp(-sigma_t * distance(ray.org, isect_->position));
                trans_pdf = exp(-sigma_t * distance(ray.org, isect_->position));
                ray.org = isect_->position + ray.dir * get_intersection_epsilon(scene);
            }

        } else if (isect_) {
            ray.org = isect_->position + ray.dir * get_intersection_epsilon(scene);
        } else {
            break;
        }

        current_path_throughput *= (transmittance / trans_pdf);
        multi_trans_pdf *= trans_pdf;

        if (!scatter && isect_ && is_light(scene.shapes[isect_->shape_id])) {
            // Reach an emissive surface
            if (never_scatter) {
                radiance += current_path_throughput * emission(*isect_, -ray.dir, scene);
            } else {
                int light_id = get_area_light_id(scene.shapes[isect_->shape_id]);
                const Light &light = scene.lights[light_id];
                Vector3 light_point_position = isect_->position;
                Vector3 light_point_normal = isect_->geometric_normal;
                PointAndNormal light_point{light_point_position, light_point_normal};

                Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, light_point, nee_p_cache, scene);

                Vector3 dir_light = normalize(light_point_position - nee_p_cache);
                Real dist = distance(nee_p_cache, light_point_position);
                Real G = max(-dot(dir_light, light_point_normal), Real(0)) / (dist * dist);
                Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);

                radiance += current_path_throughput * emission(*isect_, -ray.dir, scene) * w;
            }
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) {
            // Reach maximum bounces
            break;
        }

        if (!scatter && isect_ && isect_->material_id == -1) {
            // Index-matching surface, skip through it
            if (isect_->interior_medium_id != isect_->exterior_medium_id) {
                current_medium_id = update_medium_id(*isect_, ray, current_medium_id);
            }
            bounces++;
            continue;
        }

        // Sample next direction and update path throughput
        if (scatter) {
            never_scatter = false;

            const Medium &medium = scene.media[current_medium_id];
            Real sigma_s = get_sigma_s(medium, ray.org)[0];  // Homogeneous and monochromatic medium

            // Next event estimation
            Spectrum nee = next_event_estimation_scatter(scene, ray.org, -ray.dir, current_medium_id, bounces, rng);
            radiance += current_path_throughput * sigma_s * nee;

            // Sampling phase function
            PhaseFunction phase = get_phase_function(medium);
            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase, -ray.dir, rnd_param);
            current_path_throughput *= eval(phase, -ray.dir, *next_dir_) /
                                       pdf_sample_phase(phase, -ray.dir, *next_dir_) * sigma_s;

            // Update cache
            dir_pdf = pdf_sample_phase(phase, -ray.dir, *next_dir_);
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;

            // Update ray
            ray.dir = *next_dir_;
        } else {
            // Hit a surface -- don't need to deal with this yet
            break;
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }

    return radiance;
}

inline Spectrum next_event_estimation_surface(const Scene &scene, Vector3 p, Vector3 dir_view, int current_medium_id,
                                              const Material &mat, const PathVertex &vertex, int bounces, pcg32_state &rng) {
    // Sample a point on the light source
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 p_prime = point_on_light.position;
    Vector3 dir_light = normalize(p_prime - p);
    Real dist_light = distance(p, p_prime);
    Vector3 init_p = p;

    Spectrum T_light = make_const_spectrum(1);  // Transmittance to light
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Real p_trans_dir = 1;
    while (true) {
        Ray shadow_ray{p, dir_light,
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                           distance(p_prime, p)};
        std::optional<PathVertex> isect_ = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime);
        if (isect_) {
            next_t = distance(p, isect_->position);
        }
        if (shadow_medium_id != -1) {
            const Medium &medium = scene.media[shadow_medium_id];
            Real sigma_s = get_sigma_s(medium, p)[0];  // Homogeneous and monochromatic medium
            Real sigma_a = get_sigma_a(medium, p)[0];  // Homogeneous and monochromatic medium
            Real sigma_t = sigma_s + sigma_a;          // Monochromatic medium
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }

        if (!isect_) {
            break;
        } else {
            if (isect_->material_id >= 0) {  // Getting occluded by a non index-matching surface
                return make_zero_spectrum();
            }
            shadow_bounces++;
            if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                return make_zero_spectrum();
            }
        }

        shadow_medium_id = update_medium_id(*isect_, shadow_ray, shadow_medium_id);
        p = p + next_t * shadow_ray.dir;
    }

    if (max(T_light) > 0) {
        Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) / (dist_light * dist_light);
        Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
        Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
        Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, init_p, scene);
        Spectrum contrib = T_light * G * f * L / pdf_nee;  // cos theta seems to be included in f

        Real pdf_bsdf = pdf_sample_bsdf(mat, dir_view, dir_light, vertex, scene.texture_pool) * G * p_trans_dir;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
        return w * contrib;
    }

    return make_zero_spectrum();
}

// The fifth volumetric renderer:
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implement this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = make_const_spectrum(1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    Real dir_pdf = 0;
    Vector3 nee_p_cache{0, 0, 0};
    Real multi_trans_pdf = 1;
    bool never_scatter = true;
    bool never_bsdf = true;
    while (true) {
        bool scatter = false;
        std::optional<PathVertex> isect_ = intersect(scene, ray, ray_diff);
        Real transmittance = 1;  // Monochromatic medium
        Real trans_pdf = 1;
        if (current_medium_id != -1) {  // Not in vacuum
            const Medium &medium = scene.media[current_medium_id];
            Real sigma_s = get_sigma_s(medium, ray.org)[0];  // Homogeneous and monochromatic medium
            Real sigma_a = get_sigma_a(medium, ray.org)[0];  // Homogeneous and monochromatic medium
            Real sigma_t = sigma_s + sigma_a;                // Monochromatic medium

            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;

            if (!isect_ || t < distance(ray.org, isect_->position)) {
                scatter = true;
                transmittance = exp(-sigma_t * t);
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                ray.org += t * ray.dir;
            } else {
                transmittance = exp(-sigma_t * distance(ray.org, isect_->position));
                trans_pdf = exp(-sigma_t * distance(ray.org, isect_->position));
                ray.org = isect_->position;
            }
        } else if (isect_) {
            ray.org = isect_->position;
        } else {  // In vacuum and no intersection
            break;
        }

        current_path_throughput *= (transmittance / trans_pdf);
        multi_trans_pdf *= trans_pdf;

        if (!scatter && isect_ && is_light(scene.shapes[isect_->shape_id])) {
            // Reach an emissive surface
            if (never_scatter && never_bsdf) {
                radiance += current_path_throughput * emission(*isect_, -ray.dir, scene);
            } else {
                int light_id = get_area_light_id(scene.shapes[isect_->shape_id]);
                const Light &light = scene.lights[light_id];
                Vector3 light_point_position = isect_->position;
                Vector3 light_point_normal = isect_->geometric_normal;
                PointAndNormal light_point{light_point_position, light_point_normal};

                Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, light_point, nee_p_cache, scene);

                Vector3 dir_light = normalize(light_point_position - nee_p_cache);
                Real dist = distance(nee_p_cache, light_point_position);
                Real G = max(-dot(dir_light, light_point_normal), Real(0)) / (dist * dist);
                Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);

                radiance += current_path_throughput * emission(*isect_, -ray.dir, scene) * w;
            }
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) {
            // Reach maximum bounces
            break;
        }

        if (!scatter && isect_ && isect_->material_id == -1) {
            // Index-matching surface, skip through it
            if (isect_->interior_medium_id != isect_->exterior_medium_id) {
                current_medium_id = update_medium_id(*isect_, ray, current_medium_id);
            }
            ray.org += ray.dir * get_intersection_epsilon(scene);
            bounces++;
            continue;
        }

        // Sample next direction and update path throughput
        if (scatter) {
            never_scatter = false;

            const Medium &medium = scene.media[current_medium_id];
            Real sigma_s = get_sigma_s(medium, ray.org)[0];  // Homogeneous and monochromatic medium

            // Next event estimation
            Spectrum nee = next_event_estimation_scatter(scene, ray.org, -ray.dir, current_medium_id, bounces, rng);
            radiance += current_path_throughput * sigma_s * nee;

            // Sampling phase function
            PhaseFunction phase = get_phase_function(medium);
            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase, -ray.dir, rnd_param);
            current_path_throughput *= eval(phase, -ray.dir, *next_dir_) /
                                       pdf_sample_phase(phase, -ray.dir, *next_dir_) * sigma_s;

            // Update cache
            dir_pdf = pdf_sample_phase(phase, -ray.dir, *next_dir_);
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;

            // Update ray
            ray.dir = *next_dir_;
            ray.org += ray.dir * get_intersection_epsilon(scene);
        } else {  // Hit a non index-matching surface
            never_bsdf = false;

            const Material &mat = scene.materials[isect_->material_id];
            Vector3 dir_view = -ray.dir;

            // Next event estimation
            Spectrum nee = next_event_estimation_surface(
                scene, ray.org, dir_view, current_medium_id, mat, *isect_, bounces, rng);
            radiance += current_path_throughput * nee;

            // Sampling BSDF
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(mat,
                            dir_view,
                            *isect_,
                            scene.texture_pool,
                            bsdf_rnd_param_uv,
                            bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                break;
            }
            Vector3 dir_bsdf = bsdf_sample_->dir_out;
            Real bsdf_pdf = pdf_sample_bsdf(mat,
                                            dir_view,
                                            dir_bsdf,
                                            *isect_,
                                            scene.texture_pool);
            if (bsdf_pdf <= 0) {
                break;
            }
            Spectrum f = eval(mat,
                              dir_view,
                              dir_bsdf,
                              *isect_,
                              scene.texture_pool);
            current_path_throughput *= f / bsdf_pdf;  // cos theta seems to be included in f

            if (bsdf_sample_->eta != 0) {
                // refracted, update material
                if (isect_->interior_medium_id != isect_->exterior_medium_id) {
                    current_medium_id = update_medium_id(*isect_, ray, current_medium_id);
                }
            }

            // Update cache
            dir_pdf = bsdf_pdf;
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;

            // Update ray
            ray.dir = dir_bsdf;
            ray.org += ray.dir * get_intersection_epsilon(scene);
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }

    return radiance;
}

inline Spectrum next_event_estimation(const Scene &scene, Vector3 p, Vector3 dir_view, int current_medium_id,
                                      const Material *mat, const PathVertex *vertex, int bounces, pcg32_state &rng) {
    // Sample a point on the light source
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 p_prime = point_on_light.position;
    Vector3 dir_light = normalize(p_prime - p);
    Real dist_light = distance(p, p_prime);
    Vector3 init_p = p;

    Spectrum T_light = make_const_spectrum(1);  // Transmittance to light
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Spectrum p_trans_dir = make_const_spectrum(1);
    Spectrum p_trans_nee = make_const_spectrum(1);

    while (true) {
        Ray shadow_ray{p, dir_light,
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                           distance(p_prime, p)};
        std::optional<PathVertex> isect_ = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime);
        if (isect_) {
            next_t = distance(p, isect_->position);
        }
        if (shadow_medium_id >= 0) {
            const Medium &medium = scene.media[shadow_medium_id];
            Spectrum majorant = get_majorant(medium, shadow_ray);

            // Sample a channel for sampling
            int channel = min(2, int(next_pcg32_real<Real>(rng) * 3));
            Real accum_t = 0;
            int iteration = 0;

            // sample heterogeneous transmittance
            while (true) {
                if (majorant[channel] <= 0) {
                    break;
                }
                if (iteration >= scene.options.max_null_collisions) {
                    break;
                }

                Real u = next_pcg32_real<Real>(rng);
                Real t = -log(1 - u) / majorant[channel];
                Real dt = next_t - accum_t;
                accum_t = min(next_t, accum_t + t);
                if (t < dt) {
                    // didn't hit a surface, null-scattering event
                    p = p + t * dir_light;
                    Spectrum sigma_t = get_sigma_a(medium, p) + get_sigma_s(medium, p);
                    T_light *= exp(-majorant * t) * (majorant - sigma_t) / max(majorant);
                    p_trans_nee *= exp(-majorant * t) * majorant / max(majorant);
                    Spectrum real_prob = sigma_t / majorant;
                    p_trans_dir *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                    if (max(T_light) <= 0) {
                        break;
                    }
                } else {
                    // hit the surface
                    p = p + dt * dir_light;
                    T_light *= exp(-majorant * dt);
                    p_trans_nee *= exp(-majorant * dt);
                   
                    break;
                }
                iteration++;
            }
        }

        if (!isect_) {
            break;
        } else {
            if (isect_->material_id >= 0) {  // Getting occluded by a non index-matching surface
                return make_zero_spectrum();
            }
            shadow_bounces++;
            if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                return make_zero_spectrum();
            }

            shadow_medium_id = update_medium_id(*isect_, shadow_ray, shadow_medium_id);
            p = isect_->position;
        }
    }

    if (max(T_light) > 0) {
        T_light.x = T_light.x > 0 ? T_light.x : 0;
        T_light.y = T_light.y > 0 ? T_light.y : 0;
        T_light.z = T_light.z > 0 ? T_light.z : 0;
        if (mat) {
            if (!vertex) {
                return make_zero_spectrum();
            }

            // bsdf
            Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) / (dist_light * dist_light);
            Spectrum f = eval(*mat, dir_view, dir_light, *vertex, scene.texture_pool);
            Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
            Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, init_p, scene) * avg(p_trans_nee);
            Spectrum contrib = T_light * G * f * L / pdf_nee;  // cos theta seems to be included in f
            Real pdf_bsdf = pdf_sample_bsdf(*mat, dir_view, dir_light, *vertex, scene.texture_pool) * G * avg(p_trans_dir);
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
            return w * contrib;
        } else {
            // phase
            Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) / (dist_light * dist_light);
            PhaseFunction phase = get_phase_function(scene.media[current_medium_id]);
            Spectrum f = eval(phase, dir_view, dir_light);
            Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
            Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, init_p, scene) * avg(p_trans_nee);
            Spectrum contrib = T_light * G * f * L / pdf_nee;
            Real pdf_phase = pdf_sample_phase(phase, dir_view, dir_light) * G * avg(p_trans_dir);
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
            return w * contrib;
        }
    }

    return make_zero_spectrum();
}

// The final volumetric renderer:
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implement this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = make_const_spectrum(1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    bool never_scatter = true;
    bool never_bsdf = true;

    Real dir_pdf = 0;
    Vector3 nee_p_cache{0, 0, 0};
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    Spectrum multi_nee_pdf = make_const_spectrum(1);

    while (true) {
        bool scatter = false;
        std::optional<PathVertex> isect_ = intersect(scene, ray);
        Real t_hit = infinity<Real>();
        if (isect_) t_hit = distance(ray.org, isect_->position);
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_dir_pdf = make_const_spectrum(1);
        Spectrum trans_nee_pdf = make_const_spectrum(1);

        if (current_medium_id >= 0) {  // Not in vacuum
            const Medium &medium = scene.media[current_medium_id];
            Spectrum majorant = get_majorant(medium, ray);

            // Sample a channel for sampling
            int channel = std::clamp(int(next_pcg32_real<Real>(rng) * 3), 0, 2);
            Real accum_t = 0;
            int iteration = 0;

            // sample heterogeneous transmittance
            while (true) {
                if (majorant[channel] <= 0) {
                    break;
                }
                if (iteration >= scene.options.max_null_collisions) {
                    break;
                }

                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = t_hit - accum_t;
                accum_t = min(t_hit, accum_t + t);
                if (t < dt) {
                    // haven't reached the/a surface, sample from real/fake particle events
                    Vector3 p = ray.org + accum_t * ray.dir;
                    Spectrum sigma_t = get_sigma_a(medium, p) + get_sigma_s(medium, p);
                    Spectrum real_prob = sigma_t / majorant;
                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) {
                        // real, scatter
                        scatter = true;
                        transmittance *= exp(-majorant * t) / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * real_prob / max(majorant);
                        ray.org = ray.org + accum_t * ray.dir;
                        break;
                    } else {
                        // fake
                        transmittance *= exp(-majorant * t) * (majorant - sigma_t) / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                        trans_nee_pdf *= exp(-majorant * t) * majorant / max(majorant);
                    }
                } else {
                    // hit a surface
                    ray.org += t_hit * ray.dir;
                    transmittance *= exp(-majorant * dt);
                    trans_dir_pdf *= exp(-majorant * dt);
                    trans_nee_pdf *= exp(-majorant * dt);
                    break;
                }
                iteration++;
            }

            multi_trans_pdf *= trans_dir_pdf;
            multi_nee_pdf *= trans_nee_pdf;
            if (!scatter && !isect_) break;  // handle the case where we always hit fake particles until the last iteration
        } else if (isect_) {
            ray.org = isect_->position;
        } else {  // In vacuum and no intersection
            break;
        }

        current_path_throughput *= (transmittance / avg(trans_dir_pdf));

        if (!scatter && isect_ && is_light(scene.shapes[isect_->shape_id])) {
            // Reach an emissive surface
            if (never_scatter && never_bsdf) {
                radiance += current_path_throughput * emission(*isect_, -ray.dir, scene);
            } else {
                int light_id = get_area_light_id(scene.shapes[isect_->shape_id]);
                const Light &light = scene.lights[light_id];
                Vector3 light_point_position = isect_->position;
                Vector3 light_point_normal = isect_->geometric_normal;
                PointAndNormal light_point{light_point_position, light_point_normal};

                Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, light_point, nee_p_cache, scene) * avg(multi_nee_pdf);

                Vector3 dir_light = normalize(light_point_position - nee_p_cache);
                Real dist = distance(nee_p_cache, light_point_position);
                Real G = max(-dot(dir_light, light_point_normal), Real(0)) / (dist * dist);
                Real dir_pdf_ = dir_pdf * avg(multi_trans_pdf) * G;
                Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);

                radiance += current_path_throughput * emission(*isect_, -ray.dir, scene) * w;
            }
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) {
            // Reach maximum bounces
            break;
        }

        if (!scatter && isect_ && isect_->material_id == -1) {
            // Index-matching surface, skip through it
            if (isect_->interior_medium_id != isect_->exterior_medium_id) {
                current_medium_id = update_medium_id(*isect_, ray, current_medium_id);
            }
            Vector3 N = dot(ray.dir, isect_->geometric_normal) > 0 ? isect_->geometric_normal : -isect_->geometric_normal;
            ray.org = isect_->position + N * get_intersection_epsilon(scene);
            bounces++;
            continue;
        }

        // Sample next direction and update path throughput
        if (scatter) {
            never_scatter = false;

            const Medium &medium = scene.media[current_medium_id];
            Spectrum sigma_s = get_sigma_s(medium, ray.org);  // Homogeneous and monochromatic medium

            // Next event estimation
            Spectrum nee = next_event_estimation(scene, ray.org, -ray.dir, current_medium_id, nullptr, nullptr, bounces, rng);
            radiance += current_path_throughput * sigma_s * nee;

            // Sampling phase function
            PhaseFunction phase = get_phase_function(medium);
            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase, -ray.dir, rnd_param);
            current_path_throughput *= eval(phase, -ray.dir, *next_dir_) /
                                       pdf_sample_phase(phase, -ray.dir, *next_dir_) * sigma_s;

            // Update cache
            dir_pdf = pdf_sample_phase(phase, -ray.dir, *next_dir_);
            nee_p_cache = ray.org;
            multi_trans_pdf = make_const_spectrum(1);
            multi_nee_pdf = make_const_spectrum(1);

            // Update ray
            ray.dir = *next_dir_;
        } else {  // Hit a non index-matching surface
            never_bsdf = false;
            if (!isect_.has_value()) {
                break;
            }

            const Material &mat = scene.materials[isect_->material_id];
            Vector3 dir_view = -ray.dir;

            // Next event estimation
            Spectrum nee = next_event_estimation(
                scene, ray.org, dir_view, current_medium_id, &mat, &*isect_, bounces, rng);
            radiance += current_path_throughput * nee;

            // Sampling BSDF
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(mat,
                            dir_view,
                            *isect_,
                            scene.texture_pool,
                            bsdf_rnd_param_uv,
                            bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                break;
            }
            Vector3 dir_bsdf = bsdf_sample_->dir_out;
            Real bsdf_pdf = pdf_sample_bsdf(mat,
                                            dir_view,
                                            dir_bsdf,
                                            *isect_,
                                            scene.texture_pool);
            if (bsdf_pdf <= 0) {
                break;
            }
            Spectrum f = eval(mat,
                              dir_view,
                              dir_bsdf,
                              *isect_,
                              scene.texture_pool);
            current_path_throughput *= f / bsdf_pdf;  // cos theta seems to be included in f

            if (bsdf_sample_->eta != 0) {
                // refracted, update material
                if (isect_->interior_medium_id != isect_->exterior_medium_id) {
                    current_medium_id = update_medium_id(*isect_, ray, current_medium_id);
                }
            }

            // Update cache
            dir_pdf = bsdf_pdf;
            nee_p_cache = ray.org;
            multi_trans_pdf = make_const_spectrum(1);
            multi_nee_pdf = make_const_spectrum(1);

            // Update ray
            ray.dir = dir_bsdf;
            Vector3 N = dot(ray.dir, isect_->geometric_normal) > 0 ? isect_->geometric_normal : -isect_->geometric_normal;
            ray.org = isect_->position + N * get_intersection_epsilon(scene);
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }

    return radiance;
}