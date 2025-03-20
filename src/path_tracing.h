#pragma once

#include "scene.h"
#include "pcg.h"
#include "intersection.h"
#include <type_traits>
#include <iostream>

#define CHECK_TYPE(variable, type) (std::is_same<decltype(variable), type>::value)

///Original path tracing : 
/// Unidirectional path tracing
Spectrum path_tracing(const Scene &scene,
                      int x, int y, /* pixel coordinates */
                      pcg32_state &rng) {
    
    
    int w = scene.camera.width, h = scene.camera.height;


    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);


    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);


    // [1] HIT NOTHING OR HIT ENV MAP....
    if (!vertex_) {
        // Hit background. Account for the environment map if needed.
        if (has_envmap(scene)) {
            const Light &envmap = get_envmap(scene);
            return emission(envmap,
                            -ray.dir, // pointing outwards from light
                            ray_diff.spread,
                            PointAndNormal{}, // dummy parameter for envmap
                            scene);
        }
        return make_zero_spectrum();
    }

    // HIT SOMETHING
    PathVertex vertex = *vertex_;

    Spectrum radiance = make_zero_spectrum();
    // A path's contribution is 
    // C(v) = W(v0, v1) * G(v0, v1) * f(v0, v1, v2) * 
    //                    G(v1, v2) * f(v1, v2, v3) * 
    //                  ........
    //                  * G(v_{n-1}, v_n) * L(v_{n-1}, v_n)
    // where v is the path vertices, W is the sensor response
    // G is the geometry term, f is the BSDF, L is the emission
    //
    // "sample_primary" importance samples both W and G,
    // and we assume it always has weight 1.

    // current_path_throughput stores the ratio between
    // 1) the path contribution from v0 up to v_{i} (the BSDF f(v_{i-1}, v_i, v_{i+1}) is not included), 
    // where i is where the PathVertex "vertex" lies on, and
    // 2) the probability density for computing the path v from v0 up to v_i,
    // so that we can compute the Monte Carlo estimates C/p. 
    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    
    
    
    // eta_scale stores the scale introduced by Snell-Descartes law to the BSDF (eta^2).
    // We use the same Russian roulette strategy as Mitsuba/pbrt-v3
    // and tracking eta_scale and removing it from the
    // path contribution is crucial for many bounces of refraction.
    Real eta_scale = Real(1);

    // We hit a light immediately. 
    // This path has only two vertices and has contribution
    // C = W(v0, v1) * G(v0, v1) * L(v0, v1)
    if (is_light(scene.shapes[vertex.shape_id])) {
        radiance += current_path_throughput *
            emission(vertex, -ray.dir, scene);
    }

    // We iteratively sum up path contributions from paths with different number of vertices
    // If max_depth == -1, we rely on Russian roulette for path termination.
    int max_depth = scene.options.max_depth;
    for (int num_vertices = 3; max_depth == -1 || num_vertices <= max_depth + 1; num_vertices++) {
        // We are at v_i, and all the path contribution on and before has been accounted for.
        // Now we need to somehow generate v_{i+1} to account for paths with more vertices.
        // In path tracing, we generate two vertices:
        // 1) we sample a point on the light source (often called "Next Event Estimation")
        // 2) we randomly trace a ray from the surface point at v_i and hope we hit something.
        //
        // The first importance samples L(v_i, v_{i+1}), and the second
        // importance samples f(v_{i-1}, v_i, v_{i+1}) * G(v_i, v_{i+1})
        //
        // We then combine the two sampling strategies to estimate the contribution using weighted average.
        // Say the contribution of the first sampling is C1 (with probability density p1), 
        // and the contribution of the second sampling is C2 (with probability density p2,
        // then we compute the estimate as w1*C1/p1 + w2*C2/p2.
        //
        // Assuming the vertices for C1 is v^1, and v^2 for C2,
        // Eric Veach showed that it is a good idea setting 
        // w1 = p_1(v^1)^k / (p_1(v^1)^k + p_2(v^1)^k)
        // w2 = p_2(v^2)^k / (p_1(v^2)^k + p_2(v^2)^k),
        // where k is some scalar real number, and p_a(v^b) is the probability density of generating
        // vertices v^b using sampling method "a".
        // We will set k=2 as suggested by Eric Veach.

        // Finally, we set our "next vertex" in the loop to the v_{i+1} generated
        // by the second sampling, and update current_path_throughput using
        // our hemisphere sampling.

        // Let's implement this!
        const Material &mat = scene.materials[vertex.material_id];

        //NEE
        // First, we sample a point on the light source.
        // We do this by first picking a light source, then pick a point on it.
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, vertex.position, light_uv, shape_w, scene);

        // Next, we compute w1*C1/p1. We store C1/p1 in C1.
        Spectrum C1 = make_zero_spectrum();
        Real w1 = 0;
        // Remember "current_path_throughput" already stores all the path contribution on and before v_i.
        // So we only need to compute G(v_{i}, v_{i+1}) * f(v_{i-1}, v_{i}, v_{i+1}) * L(v_{i}, v_{i+1})
        {
            // Let's first deal with C1 = G * f * L.
            // Let's first compute G.
            Real G = 0;
            Vector3 dir_light;
            // The geometry term is different between directional light sources and
            // others. Currently we only have environment maps as directional light sources.
            if (!is_envmap(light)) {
                dir_light = normalize(point_on_light.position - vertex.position);
                // If the point on light is occluded, G is 0. So we need to test for occlusion.
                // To avoid self intersection, we need to set the tnear of the ray
                // to a small "epsilon". We set the epsilon to be a small constant times the
                // scale of the scene, which we can obtain through the get_shadow_epsilon() function.
                Ray shadow_ray{vertex.position, dir_light, 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(point_on_light.position, vertex.position)};
                if (!occluded(scene, shadow_ray)) {
                    // geometry term is cosine at v_{i+1} divided by distance squared
                    // this can be derived by the infinitesimal area of a surface projected on
                    // a unit sphere -- it's the Jacobian between the area measure and the solid angle
                    // measure.
                    G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                        distance_squared(point_on_light.position, vertex.position);
                }
            } else {
                // The direction from envmap towards the point is stored in
                // point_on_light.normal.
                dir_light = -point_on_light.normal;
                // If the point on light is occluded, G is 0. So we need to test for occlusion.
                // To avoid self intersection, we need to set the tnear of the ray
                // to a small "epsilon" which we define as c_shadow_epsilon as a global constant.
                Ray shadow_ray{vertex.position, dir_light, 
                               get_shadow_epsilon(scene),
                               infinity<Real>() /* envmaps are infinitely far away */};
                if (!occluded(scene, shadow_ray)) {
                    // We integrate envmaps using the solid angle measure,
                    // so the geometry term is 1.
                    G = 1;
                }
            }

            // Before we proceed, we first compute the probability density p1(v1)
            // The probability density for light sampling to sample our point is
            // just the probability of sampling a light times the probability of sampling a point
            Real p1 = light_pmf(scene, light_id) *
                pdf_point_on_light(light, point_on_light, vertex.position, scene);

            // We don't need to continue the computation if G is 0.
            // Also sometimes there can be some numerical issue such that we generate
            // a light path with probability zero
            if (G > 0 && p1 > 0) {
                // Let's compute f (BSDF) next.
                Vector3 dir_view = -ray.dir;
                assert(vertex.material_id >= 0);
                Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);

                // Evaluate the emission
                // We set the footprint to zero since it is not fully clear how
                // to set it in this case.
                // One way is to use a roughness based heuristics, but we have multi-layered BRDFs.
                // See "Real-time Shading with Filtered Importance Sampling" from Colbert et al.
                // for the roughness based heuristics.
                Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);

                // C1 is just a product of all of them!
                C1 = G * f * L;
            
                // Next let's compute w1

                // Remember that we want to set
                // w1 = p_1(v^1)^2 / (p_1(v^1)^2 + p_2(v^1)^2)
                // Notice that all of the probability density share the same path prefix and those cancel out.
                // Therefore we only need to account for the generation of the vertex v_{i+1}.

                // The probability density for our hemispherical sampling to sample 
                Real p2 = pdf_sample_bsdf(
                    mat, dir_view, dir_light, vertex, scene.texture_pool);
                // !!!! IMPORTANT !!!!
                // In general, p1 and p2 now live in different spaces!!
                // our BSDF API outputs a probability density in the solid angle measure
                // while our light probability density is in the area measure.
                // We need to make sure that they are in the same space.
                // This can be done by accounting for the Jacobian of the transformation
                // between the two measures.
                // In general, I recommend to transform everything to area measure 
                // (except for directional lights) since it fits to the path-space math better.
                // Converting a solid angle measure to an area measure is just a
                // multiplication of the geometry term G (let solid angle be dS, area be dA,
                // we have dA/dS = G).
                p2 *= G;

                w1 = (p1*p1) / (p1*p1 + p2*p2);
                C1 /= p1;
            }
        }
        radiance += current_path_throughput * C1 * w1;


        // Second is random shooting..
        // Let's do the hemispherical sampling next.
        Vector3 dir_view = -ray.dir;
        Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
        std::optional<BSDFSampleRecord> bsdf_sample_ =
            sample_bsdf(mat,
                        dir_view,
                        vertex,
                        scene.texture_pool,
                        bsdf_rnd_param_uv,
                        bsdf_rnd_param_w);
        if (!bsdf_sample_) {
            // BSDF sampling failed. Abort the loop.
            break;
        }
        const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
        Vector3 dir_bsdf = bsdf_sample.dir_out;
        
        // Update ray differentials & eta_scale
        if (bsdf_sample.eta == 0) {
            ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
        } else {
            ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
            eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
        }

        // Trace a ray towards bsdf_dir. Note that again we have
        // to have an "epsilon" tnear to prevent self intersection.
        Ray bsdf_ray{vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>()};
        std::optional<PathVertex> bsdf_vertex = intersect(scene, bsdf_ray);

        // To update current_path_throughput
        // we need to multiply G(v_{i}, v_{i+1}) * f(v_{i-1}, v_{i}, v_{i+1}) to it
        // and divide it with the pdf for getting v_{i+1} using hemisphere sampling.
        Real G;
        if (bsdf_vertex) {
            G = fabs(dot(dir_bsdf, bsdf_vertex->geometric_normal)) /
                distance_squared(bsdf_vertex->position, vertex.position);
        } else {
            // We hit nothing, set G to 1 to account for the environment map contribution.
            G = 1;
        }

        Spectrum f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
        Real p2 = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
        if (p2 <= 0) {
            // Numerical issue -- we generated some invalid rays.
            break;
        }

        // Remember to convert p2 to area measure!
        p2 *= G;
        // note that G cancels out in the division f/p, but we still need
        // G later for the calculation of w2.

        // Now we want to check whether dir_bsdf hit a light source, and
        // account for the light contribution (C2 & w2 & p2).
        // There are two possibilities: either we hit an emissive surface,
        // or we hit an environment map.
        // We will handle them separately.
        if (bsdf_vertex && is_light(scene.shapes[bsdf_vertex->shape_id])) {
            // G & f are already computed.
            Spectrum L = emission(*bsdf_vertex, -dir_bsdf, scene);
            Spectrum C2 = G * f * L;
            // Next let's compute p1(v2): the probability of the light source sampling
            // directly drawing the point corresponds to bsdf_dir.
            int light_id = get_area_light_id(scene.shapes[bsdf_vertex->shape_id]);
            assert(light_id >= 0);
            const Light &light = scene.lights[light_id];
            PointAndNormal light_point{bsdf_vertex->position, bsdf_vertex->geometric_normal};
            Real p1 = light_pmf(scene, light_id) *
                pdf_point_on_light(light, light_point, vertex.position, scene);
            Real w2 = (p2*p2) / (p1*p1 + p2*p2);

            C2 /= p2;
            radiance += current_path_throughput * C2;
        } else if (!bsdf_vertex && has_envmap(scene)) {
            // G & f are already computed.
            const Light &light = get_envmap(scene);
            Spectrum L = emission(light,
                                  -dir_bsdf, // pointing outwards from light
                                  ray_diff.spread,
                                  PointAndNormal{}, // dummy parameter for envmap
                                  scene);
            Spectrum C2 = G * f * L;
            // Next let's compute p1(v2): the probability of the light source sampling
            // directly drawing the direction bsdf_dir.
            PointAndNormal light_point{Vector3{0, 0, 0}, -dir_bsdf}; // pointing outwards from light
            Real p1 = light_pmf(scene, scene.envmap_light_id) *
                      pdf_point_on_light(light, light_point, vertex.position, scene);
            Real w2 = (p2*p2) / (p1*p1 + p2*p2);

            C2 /= p2;
            radiance += current_path_throughput * C2*w2;
        }

        if (!bsdf_vertex) {
            // Hit nothing -- can't continue tracing.
            break;
        }

        // Update rays/intersection/current_path_throughput/current_pdf
        // Russian roulette heuristics
        Real rr_prob = 1;
        if (num_vertices - 1 >= scene.options.rr_depth) {
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            }
        }

        ray = bsdf_ray;
        vertex = *bsdf_vertex;
        current_path_throughput = current_path_throughput * (G * f) / (p2 * rr_prob);
    }
    return radiance;
}




/// Unidirectional path tracing
GraidentPTRadiance grad_path_tracing(const Scene &scene,
                                int x, int y, /* pixel coordinates */
                                pcg32_state &rng) {
    
    // Image dimensions property & random sampling for shooting primary & shift rays :
    int w = scene.camera.width, h = scene.camera.height;
    double rng_x = next_pcg32_real<Real>(rng);
    double rng_y = next_pcg32_real<Real>(rng);
    Vector2 screen_pos((x + rng_x) / w,
                       (y + rng_y) / h);


    // Sampling a ray to be shot from the camera into the scene....
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

    // final struct to hold all information...
    GraidentPTRadiance netPath;
    // 1] Original Base Path could not be sampled :
    // We return no contribution **(mimics the generatePath = false condition in small_gdpt)**
    if (!vertex_) {
        
        // Empty struct
        return netPath;
    }



    // Sampling adjacent rays for shiftPaths & computing their first hits (all use the same random numbers) :
    // x-1, y 
    Ray ray_x0 = sample_primary(scene.camera, Vector2(((x-1) + rng_x) / w,
                                                       (y + rng_y) / h));

    std::optional<PathVertex> vertex_x0_ = intersect(scene, ray_x0, ray_diff);
    
    //x, y+1
    Ray ray_y0 = sample_primary(scene.camera, Vector2( (x + rng_x) / w,
                                                       ((y+1) + rng_y) / h));
    std::optional<PathVertex> vertex_y0_ = intersect(scene, ray_y0, ray_diff);

    // x+1, y 
    Ray ray_x1 = sample_primary(scene.camera, Vector2(((x+1) + rng_x) / w,
                                                       (y + rng_y) / h));
    std::optional<PathVertex> vertex_x1_ = intersect(scene, ray_x1, ray_diff);

    // x, y-1 
    Ray ray_y1 = sample_primary(scene.camera, Vector2( (x + rng_x) / w,
                                                       ((y-1) + rng_y) / h));
    std::optional<PathVertex> vertex_y1_ = intersect(scene, ray_y1, ray_diff);
    

    //Computing jacobians for each of the shift paths (intially there is no streching)
    double jacob_x0 = 1.0, jacob_x1 = 1.0, jacob_y0 = 1.0, jacob_y1 = 1.0;

    // Extra flag variables to indicate whether the shift path is valid to evaluate & whether it has merged or not.
    bool check_path_x0 = true, check_path_x1 = true, check_path_y0 = true, check_path_y1 = true;
    bool merge_flag_x0 = false, merge_flag_x1 = false, merge_flag_y0 = false, merge_flag_y1 = false;

    

    // Merging GDPT with LaJolla path tracing logic : 
    
   

    // HIT SOMETHING (vertex holds all important information such as location, normal, object id)
    PathVertex vertex = *vertex_;


    // If shift paths don't hit anything then do not evaluate the specific shiftPath
    if(!vertex_x0_) check_path_x0 = false;
    if(!vertex_x1_) check_path_x1 = false;
    if(!vertex_y0_) check_path_y0 = false;
    if(!vertex_y1_) check_path_y1 = false;
 
    PathVertex vertex_x0 = *vertex_x0_, vertex_x1 = *vertex_x1_, vertex_y0 = *vertex_y0_, vertex_y1 = *vertex_y1_;
    
    
    // Check if the shift & base path hit same material surfaces
    if (check_path_x0)
        if(vertex_x0.material_id != vertex.material_id) check_path_x0 = false;
        
    if (check_path_x1)
        if(vertex_x1.material_id != vertex.material_id) check_path_x1 = false;
        
    if (check_path_y0)
        if(vertex_y0.material_id != vertex.material_id) check_path_y0 = false;
        
    if (check_path_y1)
        if(vertex_y1.material_id != vertex.material_id) check_path_y1 = false;
        
    
    // Path contribution components base & all shift paths
    Spectrum contrib   = fromRGB(Vector3{1, 1, 1});
    Spectrum contribX0 = fromRGB(Vector3{1, 1, 1});
    Spectrum contribY0 = fromRGB(Vector3{1, 1, 1});
    Spectrum contribX1 = fromRGB(Vector3{1, 1, 1});
    Spectrum contribY1 = fromRGB(Vector3{1, 1, 1});

    // Probablity components base & all shift paths
    double prob = 1.0, prob_x0 = 1.0, prob_x1 = 1.0, prob_y0 = 1.0, prob_y1 = 1.0;
    
    // If path possible then update these weights...
    double wX0 = 1, wY0 = 1;
    double wX1 = 1, wY1 = 1; 
    
    // Original LaJolla's approach of evaluating a Path [IGNORING NEE for now]
    // A path's contribution is 
    // C(v) = W(v0, v1) * G(v0, v1) * f(v0, v1, v2) * 
    //                    G(v1, v2) * f(v1, v2, v3) * 
    //                  ........
    //                  * G(v_{n-1}, v_n) * L(v_{n-1}, v_n)
    // where v is the path vertices, W is the sensor response
    // G is the geometry term, f is the BSDF, L is the emission
    //
    // "sample_primary" importance samples both W and G,
    // and we assume it always has weight 1.

    // current_path_throughput stores the ratio between
    // 1) the path contribution from v0 up to v_{i} (the BSDF f(v_{i-1}, v_i, v_{i+1}) is not included), 
    // where i is where the PathVertex "vertex" lies on, and
    // 2) the probability density for computing the path v from v0 up to v_i,
    // so that we can compute the Monte Carlo estimates C/p. 

    // Note we do not modify the original code componets
    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    
    // eta_scale stores the scale introduced by Snell-Descartes law to the BSDF (eta^2).
    // We use the same Russian roulette strategy as Mitsuba/pbrt-v3
    // and tracking eta_scale and removing it from the
    // path contribution is crucial for many bounces of refraction.
    Real eta_scale = Real(1);

    // We hit a light immediately. 
    // This path has only two vertices and has contribution
    // C = W(v0, v1) * G(v0, v1) * L(v0, v1)
    if (is_light(scene.shapes[vertex.shape_id])) {
        netPath.radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
        contrib = emission(vertex, -ray.dir, scene);
    }

    // Check if alternate paths hit light or not ??
    if (check_path_x0)
        if (is_light(scene.shapes[vertex_x0.shape_id])) contribX0 = emission(vertex_x0, -ray_x0.dir, scene);
        
    if (check_path_x1)
        if (is_light(scene.shapes[vertex_x1.shape_id])) contribX1 = emission(vertex_x1, -ray_x1.dir, scene);
        
    
    if (check_path_y0)
        if (is_light(scene.shapes[vertex_y0.shape_id])) contribY0 = emission(vertex_y0, -ray_y0.dir, scene);
        
    
    if (check_path_y1)
        if (is_light(scene.shapes[vertex_y1.shape_id])) contribY1 = emission(vertex_y1, -ray_y1.dir, scene);
        
    

    // We iteratively sum up path contributions from paths with different number of vertices
    // If max_depth == -1, we rely on Russian roulette for path termination.
    int max_depth = scene.options.max_depth;
    for (int num_vertices = 3; max_depth == -1 || num_vertices <= max_depth + 1; num_vertices++) {
        // We are at v_i, and all the path contribution on and before has been accounted for.
        // Now we need to somehow generate v_{i+1} to account for paths with more vertices.
        
        // In path tracing, we generate two vertices:
        // [NOT IMPLEMENTED] : 1) we sample a point on the light source (often called "Next Event Estimation")
        // [IMPLEMENTED]     : 2) we randomly trace a ray from the surface point at v_i and hope we hit something.
        //

        // Finally, we set our "next vertex" in the loop to the v_{i+1} generated
        // by the second sampling, and update current_path_throughput using
        // our hemisphere sampling.

        // Let's implement this!
        const Material &mat = scene.materials[vertex.material_id];

        // Second is random shooting..
        // Let's do the hemispherical sampling next.

        // Hemispherical sampling for base path.
        Vector3 dir_view = -ray.dir;
        Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
        std::optional<BSDFSampleRecord> bsdf_sample_ = sample_bsdf(mat,
                                                                   dir_view,
                                                                   vertex,
                                                                   scene.texture_pool,
                                                                   bsdf_rnd_param_uv,
                                                                   bsdf_rnd_param_w);
        // Error in basePath end tracing....
        if (!bsdf_sample_) {
            // BSDF sampling failed. Abort the loop. Since base path is not possible
            return GraidentPTRadiance{};
        }
        const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
        Vector3 dir_bsdf = bsdf_sample.dir_out;
        
        // Update ray differentials & eta_scale (ssumption to take the same ray diff for next steps...)
        if (bsdf_sample.eta == 0) {
            ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
        } else {
            ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
            eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
        }

        // Trace a ray towards bsdf_dir. Note that again we have
        // to have an "epsilon" tnear to prevent self intersection.
        Ray bsdf_ray{vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>()};
        // next vertex...
        std::optional<PathVertex> bsdf_vertex = intersect(scene, bsdf_ray); 
        std::optional<PathVertex> tmp_x0_vertex = intersect(scene, Ray{vertex_x0.geometric_normal, normalize(vertex.geometric_normal - vertex_x0.geometric_normal)});
        std::optional<PathVertex> tmp_x1_vertex = intersect(scene, Ray{vertex_x1.geometric_normal, normalize(vertex.geometric_normal - vertex_x1.geometric_normal)});
        std::optional<PathVertex> tmp_y0_vertex = intersect(scene, Ray{vertex_y0.geometric_normal, normalize(vertex.geometric_normal - vertex_y0.geometric_normal)});
        std::optional<PathVertex> tmp_y1_vertex = intersect(scene, Ray{vertex_y1.geometric_normal, normalize(vertex.geometric_normal - vertex_y1.geometric_normal)});

        // Check should we merge the shift & base paths or not along with jacobians :
        if(check_path_x0)
        {
            // Both base path and shift path report lambertian surfaces....
            const Material &shift_x0_mat = scene.materials[vertex_x0.material_id];
            
            if (CHECK_TYPE(shift_x0_mat, Lambertian) and CHECK_TYPE(mat, Lambertian))
            {
                if((!tmp_x0_vertex) || (tmp_x0_vertex->shape_id != bsdf_vertex->shape_id))
                {
                    check_path_x0 = false;
                    contribX0 = make_zero_spectrum();
                    jacob_x0 = 1.0;
                }
                else
                if(!merge_flag_x0)
                {
                    merge_flag_x0 = true;
                    // If we merge it into basePath no need to evaluate this component further
                    check_path_x0 = false;
                    Vector3 baseP0 = vertex.position;
                    Vector3 p1 = bsdf_vertex->position;
                    Vector3 baseN0 = vertex.geometric_normal;
                    Vector3 n1 = bsdf_vertex->geometric_normal;

                    Vector3 baseDir = p1 - baseP0;
                    double baseDist2 = dot(baseDir,baseDir);
                    baseDir = baseDir * (1.0 / sqrt(baseDist2));
                    double baseGeom = fabs(dot(baseDir,n1)) * fabs(dot(baseDir,baseN0)) / baseDist2;
                    Vector3 shiftDir = p1 - vertex_x0.position;
                    double shiftDist2 = dot(shiftDir,shiftDir);     
                    shiftDir = shiftDir * (1.0 / sqrt(shiftDist2));
                    double shiftGeom = fabs(dot(shiftDir,n1)) * fabs(dot(shiftDir,vertex_x0.geometric_normal)) / shiftDist2;
                    jacob_x0*= (shiftGeom / baseGeom);
                }
            }

            else if(vertex_x0.material_id != vertex.material_id)
            { 
                check_path_x0 = false;
                contribX0 = make_zero_spectrum();
                jacob_x0 = 1.0;
            }
            
        }

        if(check_path_x1)
        {
            // Both base path and shift path report lambertian surfaces....
            const Material &shift_x1_mat = scene.materials[vertex_x1.material_id];
            if (CHECK_TYPE(shift_x1_mat, Lambertian) and CHECK_TYPE(mat, Lambertian))
            {
                
                if((!tmp_x1_vertex) || (tmp_x1_vertex->shape_id != bsdf_vertex->shape_id))
                {
                    check_path_x1 = false;
                    contribX1 = make_zero_spectrum();
                    jacob_x1 = 1.0;
                }
                else
                if(!merge_flag_x1)
                {
                    merge_flag_x1 = true;
                    check_path_x1 = false;
                    Vector3 baseP0 = vertex.position;
                    Vector3 p1 = bsdf_vertex->position;
                    Vector3 baseN0 = vertex.geometric_normal;
                    Vector3 n1 = bsdf_vertex->geometric_normal;

                    Vector3 baseDir = p1 - baseP0;
                    double baseDist2 = dot(baseDir,baseDir);
                    baseDir = baseDir * (1.0 / sqrt(baseDist2));
                    double baseGeom = fabs(dot(baseDir,n1)) * fabs(dot(baseDir,baseN0)) / baseDist2;
                    Vector3 shiftDir = p1 - vertex_x1.position;
                    double shiftDist2 = dot(shiftDir,shiftDir);     
                    shiftDir = shiftDir * (1.0 / sqrt(shiftDist2));
                    double shiftGeom = fabs(dot(shiftDir,n1)) * fabs(dot(shiftDir,vertex_x1.geometric_normal)) / shiftDist2;
                    jacob_x1*= (shiftGeom / baseGeom);
                }
            }
            else if(vertex_x1.material_id != vertex.material_id)
            { 
                check_path_x1 = false;
                contribX1 = make_zero_spectrum();
                jacob_x1 = 1.0;
            }
        }

        if(check_path_y0)
        {
            // Both base path and shift path report lambertian surfaces....
            const Material &shift_y0_mat = scene.materials[vertex_y0.material_id];
            if (CHECK_TYPE(shift_y0_mat, Lambertian) and CHECK_TYPE(mat, Lambertian))
            {
                
                if((!tmp_y0_vertex) || (tmp_y0_vertex->shape_id != bsdf_vertex->shape_id))
                {
                    check_path_y0 = false;
                    contribY0 = make_zero_spectrum();
                    jacob_y0 = 1.0;
                }
                else
                if(!merge_flag_y0)
                {
                    merge_flag_y0 = true;
                    check_path_y0 = false;
                    Vector3 baseP0 = vertex.position;
                    Vector3 p1 = bsdf_vertex->position;
                    Vector3 baseN0 = vertex.geometric_normal;
                    Vector3 n1 = bsdf_vertex->geometric_normal;

                    Vector3 baseDir = p1 - baseP0;
                    double baseDist2 = dot(baseDir,baseDir);
                    baseDir = baseDir * (1.0 / sqrt(baseDist2));
                    double baseGeom = fabs(dot(baseDir,n1)) * fabs(dot(baseDir,baseN0)) / baseDist2;
                    Vector3 shiftDir = p1 - vertex_x1.position;
                    double shiftDist2 = dot(shiftDir,shiftDir);     
                    shiftDir = shiftDir * (1.0 / sqrt(shiftDist2));
                    double shiftGeom = fabs(dot(shiftDir,n1)) * fabs(dot(shiftDir,vertex_x1.geometric_normal)) / shiftDist2;
                    jacob_y0*= (shiftGeom / baseGeom);
                }
            }
            else if(vertex_y0.material_id != vertex.material_id)
            { 
                check_path_y0 = false;
                contribY0 = make_zero_spectrum();
                jacob_y0 = 1.0;
            }
        }

        if(check_path_y1)
        {
            // Both base path and shift path report lambertian surfaces....
            const Material &shift_y1_mat = scene.materials[vertex_y1.material_id];
            if (CHECK_TYPE(shift_y1_mat, Lambertian) and CHECK_TYPE(mat, Lambertian))
            {
                
                if((!tmp_y1_vertex) || (tmp_y1_vertex->shape_id != bsdf_vertex->shape_id))
                {
                    check_path_y1 = false;
                    contribY1 = make_zero_spectrum();
                    jacob_y1 = 1.0;
                }
                else
                if(!merge_flag_y1)
                {
                    merge_flag_y1 = true;
                    check_path_y1 = false;
                    Vector3 baseP0 = vertex.position;
                    Vector3 p1 = bsdf_vertex->position;
                    Vector3 baseN0 = vertex.geometric_normal;
                    Vector3 n1 = bsdf_vertex->geometric_normal;

                    Vector3 baseDir = p1 - baseP0;
                    double baseDist2 = dot(baseDir,baseDir);
                    baseDir = baseDir * (1.0 / sqrt(baseDist2));
                    double baseGeom = fabs(dot(baseDir,n1)) * fabs(dot(baseDir,baseN0)) / baseDist2;
                    Vector3 shiftDir = p1 - vertex_x1.position;
                    double shiftDist2 = dot(shiftDir,shiftDir);     
                    shiftDir = shiftDir * (1.0 / sqrt(shiftDist2));
                    double shiftGeom = fabs(dot(shiftDir,n1)) * fabs(dot(shiftDir,vertex_x1.geometric_normal)) / shiftDist2;
                    jacob_y1*= (shiftGeom / baseGeom);
                }
            }
            else if(vertex_y1.material_id != vertex.material_id)
            { 
                check_path_y1 = false;
                contribY1 = make_zero_spectrum();
                jacob_y1 = 1.0;
            }
        }


        // To update current_path_throughput
        // we need to multiply G(v_{i}, v_{i+1}) * f(v_{i-1}, v_{i}, v_{i+1}) to it
        // and divide it with the pdf for getting v_{i+1} using hemisphere sampling.
        Real G;
        if (bsdf_vertex) {
            G = fabs(dot(dir_bsdf, bsdf_vertex->geometric_normal)) /
                distance_squared(bsdf_vertex->position, vertex.position);
        } else {
            // We hit nothing, set G to 1 to account for the environment map contribution.
            G = 1;
        }

        //Base path eval...
        Spectrum f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
        Real p2    = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);

       
        if (p2 <= 0) {
            // Numerical issue -- we generated some invalid rays.
            break;
        }

        // Remember to convert p2 to area measure!
        p2 *= G;
         
        // Accumulating results for base path...
        contrib = contrib*f*G;
        prob*=p2;


        // Shifted paths :
        if((check_path_x0) &&  (merge_flag_x0))
        {
            contribX0 = contribX0*f*G;
            prob_x0*=p2;
        }
        else
        if(check_path_x0)
        {
            // sample bsdf based on shifted paths
            std::optional<BSDFSampleRecord> bsdf_sample_x0_ = sample_bsdf(scene.materials[vertex_x0.material_id],
                                                                           -ray_x0.dir,
                                                                           vertex_x0,
                                                                           scene.texture_pool,
                                                                           bsdf_rnd_param_uv,
                                                                           bsdf_rnd_param_w);
            if(!bsdf_sample_x0_)
            {
                check_path_x0 = false;
                contribX0 = make_zero_spectrum();
                jacob_x0 = 1.0;
            }
            else
            {
            
                // computing shift pdf
                Real p2_x0 = pdf_sample_bsdf(scene.materials[vertex_x0.material_id], 
                                            -ray_x0.dir, 
                                             bsdf_sample_x0_->dir_out, 
                                             vertex_x0, 
                                             scene.texture_pool);
                if (p2_x0 <= 0.0)
                {
                    check_path_x0 = false;
                    contribX0 = make_zero_spectrum();
                    jacob_x0 = 1.0;
                    prob_x0 = 1.0;
                }
                else
                {
                    jacob_x0 *= p2/p2_x0;
                }
                Ray bsdf_ray_x0{vertex_x0.position, bsdf_sample_x0_->dir_out, get_intersection_epsilon(scene), infinity<Real>()};
                ray_x0 = bsdf_ray_x0;
            }
        }

        if((check_path_x1) &&  (merge_flag_x1))
        {
            contribX1 = contribX1*f*G;
            prob_x1*=p2;
        }
        else
        if(check_path_x1)
        {
            // sample bsdf based on shifted paths
            std::optional<BSDFSampleRecord> bsdf_sample_x1_ = sample_bsdf(scene.materials[vertex_x1.material_id],
                                                                           -ray_x1.dir,
                                                                           vertex_x1,
                                                                           scene.texture_pool,
                                                                           bsdf_rnd_param_uv,
                                                                           bsdf_rnd_param_w);
            if(!bsdf_sample_x1_)
            {
                check_path_x1 = false;
                contribX1 = make_zero_spectrum();
                jacob_x1 = 1.0;

            }
            else
            {
            
                // computing shift pdf
                Real p2_x1 = pdf_sample_bsdf(scene.materials[vertex_x1.material_id], 
                                            -ray_x1.dir, 
                                             bsdf_sample_x1_->dir_out, 
                                             vertex_x1, 
                                             scene.texture_pool);
                if (p2_x1 <= 0.0)
                {
                    check_path_x1 = false;
                    contribX1 = make_zero_spectrum();
                    jacob_x1 = 1.0;
                }
                else
                {
                    jacob_x1 *= p2/p2_x1;
                }
                Ray bsdf_ray_x1{vertex_x0.position, bsdf_sample_x1_->dir_out, get_intersection_epsilon(scene), infinity<Real>()};
                ray_x1 = bsdf_ray_x1;
            }

        }

        if((check_path_y0) &&  (merge_flag_y0))
        {
            contribY0 = contribY0*f*G;
            prob_y0*=p2;
        }
        else
        if(check_path_y0)
        {
                  // sample bsdf based on shifted paths
            std::optional<BSDFSampleRecord> bsdf_sample_y0_ = sample_bsdf(scene.materials[vertex_y0.material_id],
                                                                           -ray_y0.dir,
                                                                           vertex_y0,
                                                                           scene.texture_pool,
                                                                           bsdf_rnd_param_uv,
                                                                           bsdf_rnd_param_w);
            if(!bsdf_sample_y0_)
            {
                check_path_y0 = false;
                contribY0 = make_zero_spectrum();
                jacob_y0 = 1.0;

            }
            else
            {
            
                // computing shift pdf
                Real p2_y0 = pdf_sample_bsdf(scene.materials[vertex_y0.material_id], 
                                            -ray_y0.dir, 
                                             bsdf_sample_y0_->dir_out, 
                                             vertex_y0, 
                                             scene.texture_pool);
                if (p2_y0 <= 0.0)
                {
                    check_path_y0 = false;
                    contribY0 = make_zero_spectrum();
                    jacob_y0 = 1.0;
                }
                else
                {
                    jacob_y0 *= p2/p2_y0;
                }
                Ray bsdf_ray_y0{vertex_x0.position, bsdf_sample_y0_->dir_out, get_intersection_epsilon(scene), infinity<Real>()};
                ray_y0 = bsdf_ray_y0;
            }

        }

        if((check_path_y1) &&  (merge_flag_y1))
        {
            contribY1 = contribY1*f*G;
            prob_y1*=p2;
        }
        else
        if(check_path_y1)
        {
                  // sample bsdf based on shifted paths
            std::optional<BSDFSampleRecord> bsdf_sample_y1_ = sample_bsdf(scene.materials[vertex_y1.material_id],
                                                                           -ray_y1.dir,
                                                                           vertex_y1,
                                                                           scene.texture_pool,
                                                                           bsdf_rnd_param_uv,
                                                                           bsdf_rnd_param_w);
            if(!bsdf_sample_y1_)
            {
                check_path_y1 = false;
                contribY1 = make_zero_spectrum();
                jacob_y1 = 1.0;

            }
            else
            {
            
                // computing shift pdf
                Real p2_y1 = pdf_sample_bsdf(scene.materials[vertex_y1.material_id], 
                                            -ray_y1.dir, 
                                             bsdf_sample_y1_->dir_out, 
                                             vertex_y1, 
                                             scene.texture_pool);
                if (p2_y1 <= 0.0)
                {
                    check_path_y1 = false;
                    contribY1 = make_zero_spectrum();
                    jacob_y1 = 1.0;
                }
                else
                {
                    jacob_y1 *= p2/p2_y1;
                }
            }
            Ray bsdf_ray_y1{vertex_x0.position, bsdf_sample_y1_->dir_out, get_intersection_epsilon(scene), infinity<Real>()};
            ray_y1 = bsdf_ray_y1;

        }



        // note that G cancels out in the division f/p, but we still need
        // G later for the calculation of w2.

        // Now we want to check whether dir_bsdf hit a light source, and
        // account for the light contribution (C2 & w2 & p2).
        // There are two possibilities: either we hit an emissive surface,
        // or we hit an environment map.
        // We will handle them separately.
        if (bsdf_vertex && is_light(scene.shapes[bsdf_vertex->shape_id])) {
            // G & f are already computed.
            Spectrum L = emission(*bsdf_vertex, -dir_bsdf, scene);
            Spectrum C2 = G * f * L;
            contrib*=L;
            C2 /= p2;
            
            netPath.radiance += current_path_throughput * C2;

        }

        if (!bsdf_vertex) {
            // Hit nothing -- can't continue tracing.
            break;
        } else {
            // This is a valid vetex on the path so should do shift mapping for this vertex...
            // shift vertex and add its contribution to rdx0, rdy0.. etc..
        }

        // Update rays/intersection/current_path_throughput/current_pdf
        // Russian roulette heuristics
        Real rr_prob = 1;
        if (num_vertices - 1 >= scene.options.rr_depth) {
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            }
        }

        ray = bsdf_ray;
        vertex = *bsdf_vertex;
        current_path_throughput = current_path_throughput * (G * f) / (p2 * rr_prob);
        

        // Final updates :
        vertex_x0 = *tmp_x0_vertex;
        vertex_x1 = *tmp_x1_vertex;
        vertex_y0 = *tmp_y0_vertex;
        vertex_y1 = *tmp_y1_vertex;


    }

    // Final computations : 
    netPath.contrib = contrib;
    netPath.prob = prob;

    if(check_path_x0)
    {
        netPath.contribX0 = contribX0*jacob_x0;
        wX0 = prob / (prob + prob_x0*jacob_x0);
        netPath.wX0 = wX0;
    }

    if(check_path_x1)
    {
        netPath.contribX1 = contribX1*jacob_x1;
        wX1 = prob / (prob + prob_x1*jacob_x1);
        netPath.wX1 = wX1;
    }

    if(check_path_y0)
    {
        netPath.contribY0 = contribY0*jacob_y0;
        wY0 = prob / (prob + prob_y0*jacob_y0);
        netPath.wY0 = wY0;
    }

    if(check_path_y1)
    {
        netPath.contribY1 = contribY1*jacob_y1;
        wY1 = prob / (prob + prob_y1*jacob_y1);
        netPath.wY1 = wY1;
    }


    // This should return not just radiance but rdx0, rdy0 etc...
    return netPath;
}

