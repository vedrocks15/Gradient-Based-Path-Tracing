#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1: implement this!
    Vector3 half_vector = normalize(dir_in + dir_out);

    Real n_dot_in  = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real half_out  = dot(half_vector, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);

    Spectrum base_color_comp = eval(bsdf.base_color, 
                               vertex.uv, 
                               vertex.uv_screen_size, 
                               texture_pool);

    Real roughness = eval(bsdf.roughness, 
                          vertex.uv, 
                          vertex.uv_screen_size, 
                          texture_pool);

    Real anisotropic = eval(bsdf.anisotropic, 
                          vertex.uv, 
                          vertex.uv_screen_size, 
                          texture_pool);


    Real const_denom = 4*fabs(n_dot_in);
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    //Fresnel....
    Spectrum f_m = base_color_comp + (make_const_spectrum(1.0) - base_color_comp)*pow(1.0 - fabs(half_out),5);
    
    
    // Anisotropic GGX distribution
    // X_proj :
    Vector3 proj_vec = to_local(frame, half_vector);
    
    Real aspect  = sqrt(1 - 0.9*anisotropic);
    Real alpha_x = max(0.0001, pow(roughness,2)/aspect); 
    Real alpha_y = max(0.0001, pow(roughness,2)*aspect); 
    Real dist_const = c_PI*alpha_x*alpha_y;
    Real D = 1 /(dist_const*pow((pow(proj_vec.x/alpha_x,2)+  pow(proj_vec.y/alpha_y,2) + pow(proj_vec.z,2)),2));

    // Smiths model...
    Vector3 local_in = to_local(frame, dir_in);
    Vector3 local_out = to_local(frame, dir_out);

    Real internal_out = (pow(local_out.x*alpha_x,2) + pow(local_out.y*alpha_y,2)) / (pow(local_out.z,2));
    Real internal_in  = (pow(local_in.x*alpha_x,2) + pow(local_in.y*alpha_y,2)) / (pow(local_in.z,2));

    Real delta_out = (sqrt(1 + internal_out)-1)/2;
    Real delta_in  = (sqrt(1 + internal_in)-1)/2;

    Real G = (1/(1 + delta_in))*(1/(1 + delta_out));
  

    return (f_m*D*G)/const_denom;
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    Vector3 half_vector = normalize(dir_in + dir_out);

    Real n_dot_in  = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real half_out  = dot(half_vector, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);


    Real const_denom = 4*fabs(n_dot_in);

    Spectrum base_color_comp = eval(bsdf.base_color, 
                               vertex.uv, 
                               vertex.uv_screen_size, 
                               texture_pool);
    
    
    Real roughness = eval(bsdf.roughness, 
                          vertex.uv, 
                          vertex.uv_screen_size, 
                          texture_pool);

    Real anisotropic = eval(bsdf.anisotropic, 
                          vertex.uv, 
                          vertex.uv_screen_size, 
                          texture_pool);

    
    // X_proj :
    Vector3 proj_vec = to_local(frame, half_vector);
    
    Real aspect  = pow((1 - 0.9*anisotropic),0.5);
    Real alpha_x = max(0.0001, pow(roughness,2)/aspect); 
    Real alpha_y = max(0.0001, pow(roughness,2)*aspect); 

    Real dist_const = c_PI*alpha_x*alpha_y;
    Real D = 1 /(dist_const*pow((pow(proj_vec.x/alpha_x,2)+  pow(proj_vec.y/alpha_y,2) + pow(proj_vec.z,2)),2));

    Vector3 local_in = to_local(frame, dir_in);
    Real internal_in  = (pow(local_in.x*alpha_x,2) + pow(local_in.y*alpha_y,2)) / (pow(local_in.z,2));
    Real delta_in  = (sqrt(1 + internal_in)-1)/2;
    Real G = (1/(1 + delta_in));



    return (G*D)/const_denom;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    // Convert the incoming direction to local coordinates
    Vector3 local_dir_in = to_local(frame, dir_in);
    
    
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, 
                          vertex.uv, 
                          vertex.uv_screen_size, 
                          texture_pool);

    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    
    
    Real aspect  = sqrt(1 - 0.9*anisotropic);
    Real alpha_x = max(0.0001, pow(roughness,2)/aspect); 
    Real alpha_y = max(0.0001, pow(roughness,2)*aspect); 

    Real alpha = pow(roughness,2);

    // Custom function... (adding alpha_x & alpha_y for adding anisotropic behaviour)
    Vector3 local_micro_normal = sample_visible_normals_custom(local_dir_in, alpha, rnd_param_uv, alpha_x, alpha_y);


    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);

        
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    
    return BSDFSampleRecord{
            reflected,
            Real(0) /* eta */, roughness /* roughness */
        };


}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
