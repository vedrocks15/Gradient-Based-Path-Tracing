#pragma once
#include <iostream>
#include <limits>


// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    
    int w = scene.camera.width, h = scene.camera.height;

    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);

    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (vertex_) 
    {
        PathVertex vertex = *vertex_;
        Spectrum sigma_val_a = get_sigma_a(scene.media[vertex.exterior_medium_id], vertex.position);
        Real t_hit  = sqrt(pow(vertex.position[0] - ray.org[0],2) + pow(vertex.position[1] - ray.org[1],2) + pow(vertex.position[2] - ray.org[2],2));
        Spectrum transmittance = exp(-1.0*sigma_val_a*t_hit);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }
        return transmittance*Le;
    }

    return make_zero_spectrum();
}


//Helper function :
Spectrum L_s1(pcg32_state &rng, const Scene &scene, Spectrum p_sample, Ray ray, int object_id)
{
    Spectrum L_s_final = make_zero_spectrum();
    for(int i = 0; i < 1; i++)
    {
        // Sampling a point on light...
        Spectrum sigma_t = get_majorant(scene.media[object_id], ray);
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};      
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);  
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light = sample_point_on_light(light,
                                                              p_sample, 
                                                              light_uv, 
                                                              shape_w, 
                                                              scene);
        // handle internal fail case.....
        Real L_s_pdf_val =  pdf_point_on_light(light,
                                               point_on_light,
                                               p_sample,
                                               scene);
        L_s_pdf_val*=light_pmf(scene, light_id);
     
        PhaseFunction phase_f = get_phase_function(scene.media[object_id]);
        Vector3 diff = point_on_light.position - p_sample;
        Real norm_diff = sqrt(pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2));
        Vector3 omega_dash = normalize(diff);
        Spectrum phase_val = eval(phase_f, ray.dir, omega_dash);
        Spectrum Le = emission(light, -omega_dash, Real(0), point_on_light, scene);  

        // Occlusion test 
        Vector3 dir_light = normalize(point_on_light.position - p_sample);
        Ray shadow_ray{p_sample, 
                        dir_light, 
                        get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                        distance(point_on_light.position, p_sample)};  
        
        // If occlusion true then return 0
        Real occ_val = !occluded(scene, shadow_ray);
        Real jacob_val = (abs(dot(omega_dash,point_on_light.normal))/pow(norm_diff,2));
        Spectrum L_s1_estimate = phase_val*Le*exp(-sigma_t*norm_diff)*jacob_val*occ_val;
        
        L_s_final+=L_s1_estimate/L_s_pdf_val;
    }
    return L_s_final/1.0;
    
}
// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {

    // Homework 2: implememt this!

    // Sampling the primary camera ray...
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);

    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    // Computing if the ray hits something
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);


    // Sampling a distance t in the homogeneous volume space
    Real u = next_pcg32_real<Real>(rng);
    Spectrum sigma_s = get_sigma_s(scene.media[0], ray.org);
    Spectrum sigma_t = get_sigma_a(scene.media[0], ray.org) + get_sigma_s(scene.media[0], ray.org);

    Real t_val = (-1*log(1 - u))/sigma_t[0];
    Real t_max = ray.tfar; //std::numeric_limits<Real>::max();
    PathVertex vertex = *vertex_;
    Spectrum radiance = make_zero_spectrum();
    // If no vertex hit, this value is pointless...
    Real t_hit  = sqrt(pow(vertex.position[0] - ray.org[0],2) +
                       pow(vertex.position[1] - ray.org[1],2) + 
                       pow(vertex.position[2] - ray.org[2],2));
    
    if(vertex_)
    {
        t_max = t_hit;
    }
    
    if(t_val < t_max)
    {

        Spectrum transmitance = exp(-sigma_t*t_val);
        Spectrum trans_pdf    = exp(-sigma_t*t_val)*sigma_t;
        Spectrum p_sample     = ray.org + t_val*ray.dir;
        Spectrum L_s1_val     = make_zero_spectrum();
        L_s1_val = L_s1(rng, scene, p_sample, ray, 0);
        radiance+= (transmitance/trans_pdf)*sigma_s*L_s1_val;

    }
    else
    {
        // Direct Lighting condition if the object was hit...
        Spectrum trans_pdf = exp(-sigma_t*t_hit);
        Spectrum transmitance = exp(-sigma_t*t_hit);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }
        radiance += (transmitance/trans_pdf)*Le;
    }


    return radiance;
    
}

int update_medium(PathVertex vertex, Ray ray, int medium)
{
    if(vertex.interior_medium_id !=  vertex.exterior_medium_id)
    {
        if ((dot(ray.dir, vertex.geometric_normal)) > 0)
        {
            medium = vertex.exterior_medium_id;
        }
        else
        {
            medium = vertex.interior_medium_id;
        }
    }
    return medium;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    
    // Homework 2: implememt this!

    // Sampling original camera ray...
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    // Initial tracing related variables....
    int current_medium = scene.camera.medium_id;
    Real t_max = ray.tfar;
    Spectrum current_path_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;
    Spectrum trans_pdf    = make_const_spectrum(1.0);
    Spectrum transmitance = make_const_spectrum(1.0);
    
    // Tracing path....
    while (true)
    {
        // each new iteration set this variable to false..
        bool scatter = false;
        // Computing if any object is hit...
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = (*vertex_);
        
        if (current_medium >= 0)
        {
            // Sampling a distance t in the homogeneous volume space
            Real u = next_pcg32_real<Real>(rng);

            Spectrum sigma_t = get_sigma_s(scene.media[current_medium], ray.org) + get_sigma_a(scene.media[current_medium], ray.org);

            Real t_val = (-1*log(1 - u))/sigma_t[0];
            
            Real t_hit  = sqrt(pow(vertex.position[0] - ray.org[0],2) +
                               pow(vertex.position[1] - ray.org[1],2) + 
                               pow(vertex.position[2] - ray.org[2],2));

            if(vertex_)
            {
                t_max = t_hit;
            } 
            else
            {
                t_max = ray.tfar;
            }


            if (t_val < t_max)
            {
                transmitance = exp(-sigma_t*t_val);
                trans_pdf    = exp(-sigma_t*t_val)*sigma_t;
                scatter = true;
                // Updating camera ray to sample more t's
                ray.org = ray.org + t_val*ray.dir;

            }
            else
            {
                transmitance = exp(-sigma_t*t_hit);
                trans_pdf    = exp(-sigma_t*t_hit);
                // Updating camera ray to sample more t's
                ray.org = vertex.position + ray.dir*get_intersection_epsilon(scene);
            }
        }
        else
        if(vertex_)
        {
            trans_pdf    = make_const_spectrum(1.0);
            transmitance = make_const_spectrum(1.0);
            ray.org = vertex.position + ray.dir*get_intersection_epsilon(scene);
        }
        else
        {
            break;
        }

        current_path_throughput *= (transmitance/trans_pdf);

        // If no scattering...(include surface emission)
        if (!scatter)
        {
            Spectrum Le = make_zero_spectrum(); 
            if (is_light(scene.shapes[vertex.shape_id])) {
                Le = emission(vertex, -ray.dir, scene);
            }    
            radiance += current_path_throughput*Le;
        }

        //Loop break condition....
        if ((bounces == (max_depth-1)) & (max_depth != -1))
        {
            break;
        }
        
        // If scattering does not happen & a object is hit
        if ((!scatter) & (vertex_.has_value()))
        {
            // Index matching happens, update the medium
            if (vertex.material_id == -1)
            {
                // Passing through index matching surface counts as one bounce...
                current_medium = update_medium(vertex, ray, current_medium);
                bounces+=1;
                continue;
            }

        }
            
        // If scatter flag is set...
        if (scatter)
        {

            //Sample phase function
            PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            std::optional<Vector3> new_dir_ = sample_phase_function(phase_f, 
                                                                    -ray.dir, 
                                                                    Vector2{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)});
            Vector3 new_dir = *new_dir_;
            current_path_throughput *=eval(phase_f, -ray.dir, new_dir)/pdf_sample_phase(phase_f, -ray.dir, new_dir)*sigma_s;
            
            // Updating ray direction
            ray.dir = new_dir;
        }
        else
        {
            break;
        }

        
        
        Real rr_prob = 1.0;
        if (bounces >= scene.options.rr_depth)
        {
            // Russian roulette : 
            //
            rr_prob = min(current_path_throughput[0],0.95);
            if (next_pcg32_real<Real>(rng) >= rr_prob)
            {
                break;
            }
            else
            {
                current_path_throughput/=rr_prob;
            }
        }  
        bounces+=1;
    }   
    return radiance;
}

// Helper function for NEE
Spectrum next_event_estimation(pcg32_state &rng,
                               Spectrum p, 
                               Spectrum og_ray_dir,
                               int current_medium,
                               const Scene &scene,
                               int bounces)
{
    // Sampling a light source followed by sampling a point on the light source
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};      
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);  
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];

    // Reference point here is p i.e. the vertex that issues the NEE
    PointAndNormal point_on_light = sample_point_on_light(light,
                                                          p, 
                                                          light_uv, 
                                                          shape_w, 
                                                          scene);
    Spectrum p_prime = point_on_light.position;

    // Net sampling pdf value for light & point on light
    Real pdf_nee=light_pmf(scene, light_id);
    pdf_nee *=  pdf_point_on_light(light,
                                   point_on_light,
                                   p,
                                   scene);

    // Holder variables...
    Spectrum t_light       = make_const_spectrum(1.0);
    int shadow_medium  = current_medium;
    int shadow_bounces = 0;
    Real p_trans_dir   = 1.0;
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; // always set to this
    int max_depth = scene.options.max_depth;
    Vector3 dir_light = normalize(p_prime - p);

    // Computing the direction from reference point to the point sampled on the light...
    Vector3 reference_p = p;

    while (true)
    {
                
        // Shooting the shadow 
        Ray shadow_ray{p, 
                       dir_light, 
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                       distance(p_prime, p)}; 

        std::optional<PathVertex> shadow_vertex_ = intersect(scene, shadow_ray, ray_diff);
        PathVertex shadow_vertex = (*shadow_vertex_);
        Spectrum sigma_t = get_majorant(scene.media[shadow_medium], shadow_ray);
        Real next_t = distance(p,p_prime);

        // If shadow ray hits something...
        if (shadow_vertex_)
        {
            next_t = distance(p, shadow_vertex.position);
        }

        if(shadow_medium >= 0)
        {
            t_light*=exp(-sigma_t*next_t);
            p_trans_dir*=exp(-sigma_t[0]*next_t);
        }

        // Clear path to sampled point on light
        if(!shadow_vertex_)
        {
            break;
        }

        // Something is blocking the shadow ray path...
        // Full blocked case (direct return)
        // Example of visible = 0
        if(shadow_vertex.material_id>=0)
        {
            return make_zero_spectrum();
        }

        // If we enter index matching surface (early termination)
        shadow_bounces+=1;
        if ((max_depth != -1) & ((bounces + shadow_bounces + 1) >= max_depth))
        {
            return make_zero_spectrum();
        }

        // Updating the shadow medium
        shadow_medium  = update_medium(shadow_vertex, shadow_ray, shadow_medium);

        // Updating the reference point till the object at which we have hit..
        p = p + next_t*shadow_ray.dir;
    }

    
    if (max(t_light) > 0)
    {
        
        PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);
        Vector3 diff = p_prime - reference_p;
        Real norm_diff = sqrt(pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2));
        Vector3 omega_dash = normalize(diff);
        Spectrum rho = eval(phase_f, og_ray_dir, omega_dash);
        Spectrum L = emission(light, -omega_dash, Real(0), point_on_light, scene); 
        
        Real G = (max(-dot(omega_dash,point_on_light.normal),Real(0))/pow(norm_diff,2));

        Spectrum contrib = (t_light*G*rho*L)/pdf_nee;
        Real pdf_phase = pdf_sample_phase(phase_f, og_ray_dir, omega_dash)*G*p_trans_dir;
        Real w = (pdf_nee*pdf_nee)/(pdf_nee*pdf_nee + pdf_phase*pdf_phase);
        return contrib*w;
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
    
    // Homework 2: implememt this!
    
    // Sampling original camera ray...
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    // Initial tracing related variables....
    int current_medium = scene.camera.medium_id;

    Real t_max = ray.tfar;
    Spectrum current_path_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;
    Real dir_pdf = 0.0;
    Vector3 nee_path_cache = {};
    Real multi_trans_pdf = 1.0;
    bool never_scatter = true;

    // Path tracing loop 
    while (true)
    {
        // each new iteration set this variable to false..
        // Computing if any object is hit...
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = (*vertex_);
        Spectrum trans_pdf    = make_const_spectrum(1.0);
        Spectrum transmitance = make_const_spectrum(1.0);

        if (current_medium >= 0)
        {
            
            // Sampling a distance t in the homogeneous volume space
            Real u = next_pcg32_real<Real>(rng);
            Spectrum sigma_t = get_majorant(scene.media[current_medium], ray);
            Real t_val = (-1*log(1 - u))/sigma_t[0];
            
            Real t_hit  = sqrt(pow(vertex.position[0] - ray.org[0],2) +
                               pow(vertex.position[1] - ray.org[1],2) + 
                               pow(vertex.position[2] - ray.org[2],2));

            // If vertex is hit update t_max...
            if(vertex_)
            {
                t_max = t_hit;
            } 
            else
            {
                t_max = ray.tfar;
            }

            // Check for scattering or hitting...
            if (t_val < t_max)
            {
                // Nothing is hit here...
                transmitance  = exp(-sigma_t*t_val);
                trans_pdf     = exp(-sigma_t*t_val)*sigma_t;
                scatter       = true;
                // Moving t steps in the direction of camera ray
                ray.org = ray.org + t_val*ray.dir;

            }
            else
            {
                // We hit something....
                transmitance = exp(-sigma_t*t_hit);
                trans_pdf    = exp(-sigma_t*t_hit);
                ray.org = vertex.position + ray.dir*get_intersection_epsilon(scene);
            }
        }
        else
        if(vertex_)
        {
            ray.org = vertex.position + ray.dir*get_intersection_epsilon(scene);
        }
        else
        {
            break;
        }

        // Computing the net contribution along one ray...
        current_path_throughput *= (transmitance/trans_pdf);
        multi_trans_pdf*=trans_pdf[0];

        // If no scattering...(include surface emission)
        if (!scatter)
        {
            if (never_scatter)
            {
                // No need to account for MIS
                Spectrum Le = make_zero_spectrum(); 
                if (is_light(scene.shapes[vertex.shape_id])) {
                    Le = emission(vertex, -ray.dir, scene);
                }    
                radiance += current_path_throughput*Le;
            }
            else
            {
                // Need to account for MIS since phase sampling has hit a light source...
                if (is_light(scene.shapes[vertex.shape_id])) {
                        
                    // Getting back NEE pdf info
                    PathVertex light_point = vertex;
                    Real light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    const Light &light = scene.lights[light_id];

                    // Net sampling pdf value for light & point on light
                    // NEE_path_cache is the most recent vertex that issued NEE
                    Real pdf_nee =  pdf_point_on_light(light,
                                                        PointAndNormal{light_point.position, light_point.geometric_normal},
                                                        nee_path_cache,
                                                        scene);
                    pdf_nee*=light_pmf(scene, light_id);

                    //Getting phase pdf information
                    Vector3 diff = light_point.position - nee_path_cache;
                    Real norm_diff = sqrt(pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2));
                    Vector3 omega_dash = normalize(diff);
                    Real G = (max(-dot(omega_dash,light_point.geometric_normal), Real(0))/pow(norm_diff,2));
                    Real dir_pdf_ = dir_pdf*multi_trans_pdf*G;
                    Real w = (dir_pdf_*dir_pdf_)/(dir_pdf_*dir_pdf_ + pdf_nee*pdf_nee);

                        
                    // Hitting the light via phase sampling
                    Spectrum Le = emission(vertex, -ray.dir, scene);    
                    radiance += current_path_throughput*Le*w;
                }
            }
                
        }

        //Loop break condition....
        if ((bounces == (max_depth-1)) & (max_depth != -1))
        {
            break;
        }
            
        // If scattering does not happen & a object is hit
        if ((!scatter) & (vertex_.has_value()))
        {
            // Index matching happens, update the medium
            if (vertex.material_id == -1)
            {
                // Passing through index matching surface counts as one bounce...
                current_medium = update_medium(vertex, ray, current_medium);
                bounces+=1;
                continue;
            }
        }
                
        // If scatter flag is set...
        if (scatter)
        {
            never_scatter = false;
                
            // Sample between phase function & NEE
            // Vertex at which we decide to do NEE
            nee_path_cache = ray.org;
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum L_nee = next_event_estimation(rng, nee_path_cache, -ray.dir,current_medium, scene, bounces);

            // Adding the result of NEE 
            radiance += current_path_throughput*sigma_s*L_nee;
        
            //Sample phase function
            PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);
            std::optional<Vector3> new_dir_ = sample_phase_function(phase_f, 
                                                                    -ray.dir, 
                                                                    Vector2{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)});
            Vector3 new_dir = *new_dir_;
            current_path_throughput *=eval(phase_f, -ray.dir, new_dir)[0]/pdf_sample_phase(phase_f, -ray.dir, new_dir)*sigma_s;
                
                
            dir_pdf = pdf_sample_phase(phase_f, -ray.dir, new_dir);
            multi_trans_pdf = 1.0;

            // Updating ray direction
            ray.dir = new_dir;
                
        }
        else
        {
            // No scattering & no hit...
            break;
        }

    
        // Russian roulette : 
        Real rr_prob = 1.0;
        if (bounces >= scene.options.rr_depth)
        {
            rr_prob = min(max(current_path_throughput),0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob)
            {
                break;
            }
            else
            {
                current_path_throughput/=rr_prob;
            }
        }  
        bounces+=1;
    }   
    return radiance;
}

Spectrum next_event_estimation_surface_case(pcg32_state &rng,
                               Spectrum p, 
                               Spectrum og_ray_dir,
                               int current_medium,
                               const Scene &scene,
                               int bounces,
                               const Material &mat,
                               PathVertex vertex)
{
    // Sampling a light source followed by sampling a point on the light source
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};      
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);  
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];

    // Reference point here is p i.e. the vertex that issues the NEE
    PointAndNormal point_on_light = sample_point_on_light(light,
                                                          p, 
                                                          light_uv, 
                                                          shape_w, 
                                                          scene);
    Spectrum p_prime = point_on_light.position;
    Vector3 dir_light = normalize(p_prime - p);

    // Computing the direction from reference point to the point sampled on the light...
    Vector3 reference_p = p;

    // Net sampling pdf value for light & point on light
    Real pdf_nee=light_pmf(scene, light_id);
    pdf_nee *=  pdf_point_on_light(light,
                                   point_on_light,
                                   p,
                                   scene);

    // Holder variables...
    Spectrum t_light       = make_const_spectrum(1.0);
    int shadow_medium  = current_medium;
    int shadow_bounces = 0;
    Real p_trans_dir   = 1.0;
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; // always set to this
    int max_depth = scene.options.max_depth;
    
    while (true)
    {
                
        // Shooting the shadow 
        Ray shadow_ray{p, 
                       dir_light, 
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                       distance(p_prime, p)}; 

        std::optional<PathVertex> shadow_vertex_ = intersect(scene, shadow_ray, ray_diff);
        PathVertex shadow_vertex = (*shadow_vertex_);
        Spectrum sigma_t = get_majorant(scene.media[shadow_medium], shadow_ray);
        Real next_t = distance(p,p_prime);

        // If shadow ray hits something...
        if (shadow_vertex_)
        {
            next_t = distance(p, shadow_vertex.position);
        }

        if(shadow_medium >= 0)
        {
            t_light*=exp(-sigma_t*next_t);
            p_trans_dir*=exp(-sigma_t[0]*next_t);
        }

        // Clear path to sampled point on light
        if(!shadow_vertex_)
        {
            break;
        }

        // Something is blocking the shadow ray path...
        // Full blocked case (direct return)
        // Example of visible = 0
        if(shadow_vertex.material_id>=0)
        {
            return make_zero_spectrum();
        }

        // If we enter index matching surface (early termination)
        shadow_bounces+=1;
        if ((max_depth != -1) & ((bounces + shadow_bounces + 1) >= max_depth))
        {
            return make_zero_spectrum();
        }

        // Updating the shadow medium
        shadow_medium  = update_medium(shadow_vertex, shadow_ray, shadow_medium);

        // Updating the reference point till the object at which we have hit..
        p = p + next_t*shadow_ray.dir;
    }

    
    if (max(t_light) > 0)
    {
        
        Vector3 diff = p_prime - reference_p;
        Real norm_diff = sqrt(pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2));
        Vector3 omega_dash = normalize(diff);
        Spectrum rho = eval(mat, og_ray_dir, omega_dash, vertex, scene.texture_pool);
        Spectrum L = emission(light, -omega_dash, Real(0), point_on_light, scene); 
        
        Real G = (max(-dot(omega_dash,point_on_light.normal),Real(0))/pow(norm_diff,2));

        Spectrum contrib = (t_light*G*rho*L)/pdf_nee;
        Real pdf_surface = pdf_sample_bsdf(mat, og_ray_dir, omega_dash, vertex, scene.texture_pool)*G*p_trans_dir;
        Real w = (pdf_nee*pdf_nee)/(pdf_nee*pdf_nee + pdf_surface*pdf_surface);
        return contrib*w;
    }
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng)  {
    
    // Homework 2: implememt this!
    
    // Sampling original camera ray...
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    // Setting up inital variables before walking along the sampled ray...
    int current_medium = scene.camera.medium_id;
    Real t_max = ray.tfar;
    Spectrum current_path_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;
    Real dir_pdf = 0.0;
    Vector3 nee_path_cache = {};
    Real multi_trans_pdf = 1.0;
    bool never_scatter = true;
    bool never_surface = true;

    // Path tracing loop 
    while (true)
    {
        // each new iteration set this variable to false..
        // Computing if any object is hit...
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = (*vertex_);

        // Variables to compute details of current medium
        Spectrum trans_pdf    = make_const_spectrum(1.0);
        Spectrum transmitance = make_const_spectrum(1.0);

        if (current_medium >= 0)
        {
            
            // Sampling a distance t in the homogeneous volume space
            Real u = next_pcg32_real<Real>(rng);
            Spectrum sigma_t = get_majorant(scene.media[current_medium], ray);
            Real t_val = (-1*log(1 - u))/sigma_t[0];
            
            Real t_hit  = sqrt(pow(vertex.position[0] - ray.org[0],2) +
                               pow(vertex.position[1] - ray.org[1],2) + 
                               pow(vertex.position[2] - ray.org[2],2));

            // If vertex is hit update t_max...
            if(vertex_)
            {
                t_max = t_hit;
            } 
            else
            {
                t_max = ray.tfar;
            }

            // Check for scattering or hitting...
            if (t_val < t_max)
            {
                // Nothing is hit here...
                transmitance  = exp(-sigma_t*t_val);
                trans_pdf     = exp(-sigma_t*t_val)*sigma_t;
                scatter       = true;
                never_scatter = false;

                // Moving t steps in the direction of camera ray
                ray.org = ray.org + t_val*ray.dir;

            }
            else
            {
                // We hit something....
                transmitance = exp(-sigma_t*t_hit);
                trans_pdf    = exp(-sigma_t*t_hit);

                // Updating the ray to t_hit vertex position
                ray.org = vertex.position;
            }
        }
        else
        if(vertex_)
        {
            // Medium id -1
            ray.org = vertex.position;
        }
        else
        {
            break;
        }

        // Computing the net contribution along one ray...
        current_path_throughput *= (transmitance/trans_pdf);
        multi_trans_pdf*=trans_pdf[0];

        // If no scattering...(include surface emission)
        if (!scatter)
        {
            // Hitting the light source directly
            if ((never_scatter) && (never_surface))
            {
                // No need to account for MIS
                Spectrum Le = make_zero_spectrum(); 
                if (is_light(scene.shapes[vertex.shape_id])) {
                    Le = emission(vertex, -ray.dir, scene);
                }    
                radiance += current_path_throughput*Le;
            }
            else
            {
                // Need to account for MIS since phase sampling has hit a light source...
                if (is_light(scene.shapes[vertex.shape_id])) {
                        
                    // Getting back NEE pdf info
                    PathVertex light_point = vertex;
                    Real light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    const Light &light = scene.lights[light_id];

                    // Net sampling pdf value for light & point on light
                    // NEE_path_cache is the most recent vertex that issued NEE
                    Real pdf_nee =  pdf_point_on_light(light,
                                                        PointAndNormal{light_point.position, light_point.geometric_normal},
                                                        nee_path_cache,
                                                        scene);
                    pdf_nee*=light_pmf(scene, light_id);

                    //Getting phase pdf information
                    Vector3 diff = light_point.position - nee_path_cache;
                    Real norm_diff = sqrt(pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2));
                    Vector3 omega_dash = normalize(diff);
                    Real G = (max(-dot(omega_dash,light_point.geometric_normal), Real(0))/pow(norm_diff,2));
                    Real dir_pdf_ = dir_pdf*multi_trans_pdf*G;
                    Real w = (dir_pdf_*dir_pdf_)/(dir_pdf_*dir_pdf_ + pdf_nee*pdf_nee);

                        
                    // Hitting the light via phase sampling or via bsdf
                    // The change is accounted in dir_pdf & multi_trans_pdf...
                    Spectrum Le = emission(vertex, -ray.dir, scene);    
                    radiance += current_path_throughput*Le*w;
                }
            }
                
        }

        //Loop break condition....
        if ((bounces == (max_depth-1)) & (max_depth != -1))
        {
            break;
        }
            
        // If scattering does not happen & a object is hit
        if ((!scatter) & (vertex_.has_value()))
        {
            // Index matching happens, update the medium
            if (vertex.material_id == -1)
            {
                // Passing through index matching surface counts as one bounce...
                current_medium = update_medium(vertex, ray, current_medium);
                bounces+=1;
                continue;
            }
        }
                
        // If scatter flag is set...
        if (scatter)
        {
                
            // Sample between phase function & NEE
            // Vertex at which we decide to do NEE
            nee_path_cache = ray.org;
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum L_nee = next_event_estimation(rng, nee_path_cache, -ray.dir,current_medium, scene, bounces);

            // Adding the result of NEE 
            radiance += current_path_throughput*sigma_s*L_nee;
        
            //Sample phase function
            PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);
            std::optional<Vector3> new_dir_ = sample_phase_function(phase_f, 
                                                                    -ray.dir, 
                                                                    Vector2{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)});
            Vector3 new_dir = *new_dir_;
            current_path_throughput *=eval(phase_f, -ray.dir, new_dir)/pdf_sample_phase(phase_f, -ray.dir, new_dir)*sigma_s;
                
            // Setting up cache variables...
            dir_pdf = pdf_sample_phase(phase_f, -ray.dir, new_dir);
            multi_trans_pdf = 1.0;

            // Updating ray direction
            ray.dir = new_dir;
            // Prevent self-intersection therefore adding some epsillon
            ray.org += ray.dir * get_intersection_epsilon(scene);
                
        }
        else
        {
            // hit a surface (Hitting a surface that would add color & direction)
            // At each surface we have to take NEE & BSDF sampling component also...
            never_surface = false;
            nee_path_cache = ray.org;
            // Surface implementation of NEE...
            Spectrum L_nee = next_event_estimation_surface_case(rng, 
                                                                nee_path_cache, 
                                                                -ray.dir, 
                                                                current_medium, 
                                                                scene, 
                                                                bounces, 
                                                                scene.materials[vertex.material_id], vertex);
            radiance += current_path_throughput * L_nee;

            // Taken from path tracing code
            //Sampling bsdf
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(scene.materials[vertex.material_id],
                            -ray.dir,
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

            Real p2 = pdf_sample_bsdf(scene.materials[vertex.material_id], 
                                      -ray.dir, 
                                      dir_bsdf, 
                                      vertex, 
                                      scene.texture_pool);
            if (p2 <= 0) {
                // Numerical issue -- we generated some invalid rays.
                break;
            }
            Spectrum f = eval(scene.materials[vertex.material_id], 
                              -ray.dir, 
                              dir_bsdf, 
                              vertex, 
                              scene.texture_pool);


            current_path_throughput *= f / p2;

            if (bsdf_sample.eta != 0) {
                    // refracted, update material
                    current_medium = update_medium(vertex, ray, current_medium);
            }

            dir_pdf = p2; 
            multi_trans_pdf = 1;  
            // Update ray
            ray.dir = dir_bsdf;
            ray.org += ray.dir * get_intersection_epsilon(scene);
        } 


    
        // Russian roulette : 
        Real rr_prob = 1.0;
        if (bounces >= scene.options.rr_depth)
        {
            rr_prob = min(max(current_path_throughput),0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob)
            {
                break;
            }
            else
            {
                current_path_throughput/=rr_prob;
            }
        }  
        bounces+=1;
    }   
    return radiance;
}


// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum next_event_estimation_final(pcg32_state &rng,
                               Spectrum p, 
                               Spectrum og_ray_dir,
                               int current_medium,
                               const Scene &scene,
                               int bounces,
                               const Material *mat,
                               PathVertex *vertex)
{
    // Sampling a light source followed by sampling a point on the light source
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};      
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);  
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];

    // Reference point here is p i.e. the vertex that issues the NEE
    PointAndNormal point_on_light = sample_point_on_light(light,
                                                          p, 
                                                          light_uv, 
                                                          shape_w, 
                                                          scene);
    Spectrum p_prime  = point_on_light.position;
    Vector3 dir_light = normalize(p_prime - p);

    // Computing the direction from reference point to the point sampled on the light...
    Vector3 reference_p = p;

    // Holder variables...
    Spectrum t_light   = make_const_spectrum(1.0);
    int shadow_medium  = current_medium;
    int shadow_bounces = 0;
    Spectrum p_trans_nee = make_const_spectrum(1.0);
    Spectrum p_trans_dir = make_const_spectrum(1.0);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; // always set to this
    int max_depth = scene.options.max_depth;
    
    while (true)
    {
                
        // Shooting the shadow 
        Ray shadow_ray{p, 
                       dir_light, 
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                       distance(p_prime, p)}; 

        std::optional<PathVertex> shadow_vertex_ = intersect(scene, shadow_ray, ray_diff);
        PathVertex shadow_vertex = (*shadow_vertex_);
        
        Real next_t = distance(p,p_prime);

        // If shadow ray hits something...
        if (shadow_vertex_)
        {
            next_t = distance(p, shadow_vertex.position);
        }

        if(shadow_medium >= 0)
        {  
            Spectrum sigma_t_m = get_majorant(scene.media[shadow_medium], shadow_ray);
            // Sampling a channel : 
            int random_channel =  min(2, int(next_pcg32_real<Real>(rng) * 3));//std::clamp(int(next_pcg32_real<Real>(rng) * 3), 0, 2);
            Real accum_t = 0;
            int iteration = 0;

            while(true)
            {
                if (sigma_t_m[random_channel] <= 0) {
                    break;
                }
                if (iteration >= scene.options.max_null_collisions){
                    break;
                }

                Real u = next_pcg32_real<Real>(rng);
                Real t_val = (-1*log(1 - u))/sigma_t_m[random_channel];
                Real dt = next_t - accum_t;
                accum_t = min(accum_t + t_val, next_t);
                if(t_val < dt)
                {
                    p = p + t_val*shadow_ray.dir;
                    Spectrum sigma_t = get_sigma_a(scene.media[shadow_medium], p) + get_sigma_s(scene.media[shadow_medium], p);
                    t_light *= exp(-sigma_t_m * t_val) * (sigma_t_m - sigma_t) / max(sigma_t_m);
                    p_trans_nee *= exp(-sigma_t_m * t_val) * sigma_t_m / max(sigma_t_m);
                    Spectrum real_prob = sigma_t / sigma_t_m;
                    p_trans_dir *= exp(-sigma_t_m * t_val) * sigma_t_m * (1 - real_prob) / max(sigma_t_m);
                    if (max(t_light) <= 0)
                    {
                        break;
                    }
                }
                else
                {
                    p = p + dt*shadow_ray.dir;
                    t_light *= exp(-sigma_t_m * dt);
                    p_trans_nee *= exp(-sigma_t_m * dt);
                    p_trans_dir *= exp(-sigma_t_m * dt);
                    break;
                }
                iteration+=1;
            }
        }
        // Clear path to sampled point on light
        if(!shadow_vertex_)
        {
            break;
        }

        // Something is blocking the shadow ray path...
        // Full blocked case (direct return)
        // Example of visible = 0
        if(shadow_vertex.material_id>=0)
        {
            return make_zero_spectrum();
        }

        // If we enter index matching surface (early termination)
        shadow_bounces+=1;
        if ((max_depth != -1) & ((bounces + shadow_bounces + 1) >= max_depth))
        {
            return make_zero_spectrum();
        }

        // Updating the shadow medium
        shadow_medium  = update_medium(shadow_vertex, shadow_ray, shadow_medium);

        // Updating the reference point till the object at which we have hit..
        p = shadow_vertex.position;
    }

    // Combining NEE of phase sampling & bsdf
    if (max(t_light) > 0)
    {
        // BSDF CASE
        if(mat)
        {
            if (!vertex)
            {
                return make_zero_spectrum();
            }
            
            // Net sampling pdf value for light & point on light
            Real pdf_nee=light_pmf(scene, light_id);
            pdf_nee *=  pdf_point_on_light(light,
                                           point_on_light,
                                           reference_p,
                                           scene);
            pdf_nee *= avg(p_trans_nee);

            // bsdf case 
            Vector3 diff = p_prime - reference_p;
            Real norm_diff = sqrt(pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2));
            Vector3 omega_dash = normalize(diff);
            Spectrum rho = eval(*mat, og_ray_dir, omega_dash, *vertex, scene.texture_pool);
            Spectrum L = emission(light, -omega_dash, Real(0), point_on_light, scene); 
            
            Real G = (max(-dot(omega_dash,point_on_light.normal),Real(0))/pow(distance(p_prime, reference_p),2));

            Spectrum contrib = (t_light*G*rho*L)/pdf_nee;
            Real pdf_bsdf = pdf_sample_bsdf(*mat, og_ray_dir, omega_dash, *vertex, scene.texture_pool)*G*avg(p_trans_dir);
            Real w = (pdf_nee*pdf_nee)/(pdf_nee*pdf_nee + pdf_bsdf*pdf_bsdf);
            return contrib*w;

        }
        else
        {
            // SCATTERING CASE
            Real pdf_nee=light_pmf(scene, light_id);
            pdf_nee *=  pdf_point_on_light(light,
                                        point_on_light,
                                        reference_p,
                                        scene);
            pdf_nee *= avg(p_trans_nee);
            
            PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);
            Vector3 diff = p_prime - reference_p;
            Real norm_diff = sqrt(pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2));
            Vector3 omega_dash = normalize(diff);
            Spectrum rho = eval(phase_f, og_ray_dir, omega_dash);
            Spectrum L = emission(light, -omega_dash, Real(0), point_on_light, scene); 
            
            Real G = (max(-dot(omega_dash,point_on_light.normal),Real(0))/pow(norm_diff,2));

            Spectrum contrib = (t_light*G*rho*L)/pdf_nee;
            Real pdf_phase = pdf_sample_phase(phase_f, og_ray_dir, omega_dash)*G*avg(p_trans_dir);
            Real w = (pdf_nee*pdf_nee)/(pdf_nee*pdf_nee + pdf_phase*pdf_phase);
            return contrib*w;
    
        }
        
        
    }
    return make_zero_spectrum();
}


Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    
    // Homework 2: implememt this!
    
    // Sampling original camera ray...
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    // Setting inital variables & caching ones too...
    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;
    Vector3 nee_path_cache = {};
    bool never_scatter = true;
    bool never_surface = true;
    Real dir_pdf = 0.0;
    Spectrum multi_trans_pdf     = make_const_spectrum(1.0);
    Spectrum mult_trans_nee_pdf  = make_const_spectrum(1.0);

    // Path tracing loop 
    while (true)
    {
        // each new iteration set this variable to false..
        // Computing if any object is hit...
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = (*vertex_);
        Spectrum transmitance  = make_const_spectrum(1.0);
        Spectrum trans_dir_pdf = make_const_spectrum(1.0);
        Spectrum trans_nee_pdf = make_const_spectrum(1.0);
        Spectrum sigma_t_m =  get_majorant(scene.media[current_medium], ray);
        
        // Check if there was a hit or not
        Real t_hit = 0.0;
        if(vertex_)
        {
            t_hit = sqrt(pow(vertex.position[0] - ray.org[0],2) +
                         pow(vertex.position[1] - ray.org[1],2) + 
                         pow(vertex.position[2] - ray.org[2],2));
        }
        else
        {
            t_hit = ray.tfar;
        }

        if (current_medium >= 0)
        {
            
            // Sampling a channel : (for chromatic case) 
            int random_channel = std::clamp(int(next_pcg32_real<Real>(rng) * 3), 0, 2);
            Real accum_t = 0;
            int iteration = 0;

            // while loop added so that we keep on moving till we hit a real particle...
            while(true)
            {
                if (sigma_t_m[random_channel] <= 0) {
                    break;
                }
                if (iteration >= scene.options.max_null_collisions){
                    break;
                }

                // Sampled distance...
                Real u = next_pcg32_real<Real>(rng);
                Real t_val = (-1*log(1 - u))/sigma_t_m[random_channel];
                // Remaining distacne to hit
                Real dt = t_hit - accum_t;

                // Total completed distance & capping this to t_hit so that dt is never < 0
                accum_t = min(accum_t + t_val, t_hit);

                if (t_val < dt)
                {
                    // Distance reached...
                    Vector3 reached_point = ray.org + accum_t*ray.dir;
                    Spectrum sigma_t = get_sigma_a(scene.media[current_medium], reached_point) + get_sigma_s(scene.media[current_medium], reached_point);
                    
                    // Fake or real particle...
                    Spectrum real_prob = sigma_t / sigma_t_m;
                    if(next_pcg32_real<Real>(rng) < real_prob[random_channel])
                    {
                        //Scatter Real particle.... (end)
                        scatter = true;
                        never_scatter = false;
                        transmitance *= exp(-sigma_t_m * t_val) / max(sigma_t_m);
                        trans_dir_pdf *= exp(-sigma_t_m * t_val) * sigma_t_m * real_prob / max(sigma_t_m);
                        //don’t need to account for trans_nee_pdf since we scatter
                        
                        // Update ray...
                        ray.org = ray.org + accum_t * ray.dir;
                        break;
                    } 
                    else
                    {
                        // Fake particle...
                        transmitance  *= exp(-sigma_t_m * t_val) * (sigma_t_m - sigma_t) / max(sigma_t_m);
                        trans_dir_pdf *= exp(-sigma_t_m * t_val) * sigma_t_m * ( 1- real_prob) / max(sigma_t_m);
                        trans_nee_pdf *= exp(-sigma_t_m * t_val) * sigma_t_m / max(sigma_t_m);
                    }
                }
                else
                {
                    // hit a real surface...
                    ray.org += t_hit * ray.dir;
                    transmitance  *= exp(-sigma_t_m*dt);
                    trans_dir_pdf *= exp(-sigma_t_m*dt);
                    trans_nee_pdf *= exp(-sigma_t_m*dt);
                    break;
                }
                
                iteration+=1;


            }

            // Free flight sampling or NEE pdf accumulation
            multi_trans_pdf *= trans_dir_pdf;
            mult_trans_nee_pdf *= trans_nee_pdf;

        }
        else
        if(vertex_)
        {
            // Direct hit update position...
            ray.org = vertex.position;
        }
        else
        {
            break;
        }

        // Computing the net contribution along one ray...
        current_path_throughput *= (transmitance/avg(trans_dir_pdf));

        // If no scattering...(include surface emission)
        if (!scatter)
        {
            if ((never_scatter) && (never_surface))
            {
                // No need to account for MIS
                Spectrum Le = make_zero_spectrum(); 
                if (is_light(scene.shapes[vertex.shape_id])) {
                    Le = emission(vertex, -ray.dir, scene);
                }    
                radiance += current_path_throughput*Le;
            }
            else
            {
                // Need to account for MIS since phase sampling has hit a light source...
                if (is_light(scene.shapes[vertex.shape_id])) {
                        
                    // Getting back NEE pdf info
                    PathVertex light_point = vertex;
                    Real light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    const Light &light = scene.lights[light_id];

                    // Net sampling pdf value for light & point on light
                    // NEE_path_cache is the most recent vertex that issued NEE
                    Real pdf_nee =  pdf_point_on_light(light,
                                                        PointAndNormal{light_point.position, light_point.geometric_normal},
                                                        nee_path_cache,
                                                        scene);

                    pdf_nee*=light_pmf(scene, light_id)*avg(mult_trans_nee_pdf);

                    //Getting phase pdf information
                    Vector3 diff = light_point.position - nee_path_cache;
                    Real norm_diff = sqrt(pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2));
                    Vector3 omega_dash = normalize(diff);
                    Real G = (max(-dot(omega_dash,light_point.geometric_normal), Real(0))/pow(norm_diff,2));
                    Real dir_pdf_ = dir_pdf*avg(multi_trans_pdf)*G;
                    Real w = (dir_pdf_*dir_pdf_)/(dir_pdf_*dir_pdf_ + pdf_nee*pdf_nee);

                        
                    // Hitting the light via phase sampling
                    Spectrum Le = emission(vertex, -ray.dir, scene);    
                    radiance += current_path_throughput*Le*w;
                }
            }
                
        }

        //Loop break condition....
        if ((bounces == (max_depth-1)) & (max_depth != -1))
        {
            break;
        }
            
        // If scattering does not happen & a object is hit
        if ((!scatter) & (vertex_.has_value()))
        {
            // Index matching happens, update the medium
            if (vertex.material_id == -1)
            {
                // Passing through index matching surface counts as one bounce...
                current_medium = update_medium(vertex, ray, current_medium);
                Vector3 final_dir = dot(ray.dir, vertex.geometric_normal) > 0 ? vertex.geometric_normal : -vertex.geometric_normal;
                ray.org = vertex.position + final_dir * get_intersection_epsilon(scene);
                bounces+=1;
                continue;
            }
        }
                
        // If scatter flag is set...
        if (scatter)
        {
            
            // Sample between phase function & NEE
            // Vertex at which we decide to do NEE
            nee_path_cache = ray.org;
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum L_nee = next_event_estimation_final(rng, 
                                                         nee_path_cache, 
                                                         -ray.dir,
                                                         current_medium, 
                                                         scene, 
                                                         bounces, 
                                                         nullptr, 
                                                         nullptr);

            // Adding the result of NEE 
            radiance += current_path_throughput*sigma_s*L_nee;
        
            //Sample phase function
            PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);
            std::optional<Vector3> new_dir_ = sample_phase_function(phase_f, 
                                                                    -ray.dir, 
                                                                    Vector2{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)});
            Vector3 new_dir = *new_dir_;
            current_path_throughput *=eval(phase_f, -ray.dir, new_dir)/pdf_sample_phase(phase_f, -ray.dir, new_dir)*sigma_s;
                
                
            dir_pdf = pdf_sample_phase(phase_f, -ray.dir, new_dir);
            multi_trans_pdf = make_const_spectrum(1);
            mult_trans_nee_pdf = make_const_spectrum(1);


            // Updating ray direction
            ray.dir = new_dir;                
        }
        else
        {
            // hit a surface
            never_surface = false;
            nee_path_cache = ray.org;

            Spectrum L_nee = next_event_estimation_final(rng, 
                                                         nee_path_cache, 
                                                         -ray.dir, 
                                                         current_medium, 
                                                         scene, 
                                                         bounces, 
                                                         &scene.materials[vertex.material_id], 
                                                         &vertex);

            radiance += current_path_throughput * L_nee;

            //Sampling bsdf
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(scene.materials[vertex.material_id],
                            -ray.dir,
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

            Real p2 = pdf_sample_bsdf(scene.materials[vertex.material_id], 
                                      -ray.dir, 
                                      dir_bsdf, 
                                      vertex, 
                                      scene.texture_pool);
            if (p2 <= 0) {
                // Numerical issue -- we generated some invalid rays.
                break;
            }
            Spectrum f = eval(scene.materials[vertex.material_id], 
                              -ray.dir, 
                              dir_bsdf, 
                              vertex, 
                              scene.texture_pool);


            current_path_throughput *= f / p2;

            if (bsdf_sample.eta != 0) {
                    // refracted, update material
                    current_medium = update_medium(vertex, ray, current_medium);
            }

            dir_pdf = p2;
            multi_trans_pdf = make_const_spectrum(1);
            mult_trans_nee_pdf = make_const_spectrum(1);

            // Update ray
            ray.dir = dir_bsdf;
            //Vector3 final_dir = dot(ray.dir, vertex.geometric_normal) > 0 ? vertex.geometric_normal : -vertex.geometric_normal;
            ray.org =  vertex.position +  ray.dir*get_intersection_epsilon(scene);
        } 

    
        // Russian roulette : 
        Real rr_prob = 1.0;
        if (bounces >= scene.options.rr_depth)
        {
            rr_prob = min(max(current_path_throughput),0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob)
            {
                break;
            }
            else
            {
                current_path_throughput/=rr_prob;
            }
        }  
        bounces+=1;
    }   
    return radiance;
}


