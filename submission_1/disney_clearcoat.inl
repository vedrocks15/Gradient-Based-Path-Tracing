#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real r_0 = pow((1.5-1.0),2) / pow((1.5+1.0),2);
    Real f_c = r_0 + (1-r_0)*pow(1 - fabs(half_out), 5);

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, 
                               vertex.uv, 
                               vertex.uv_screen_size, 
                               texture_pool);

    Real alpha_g = (1-clearcoat_gloss)*0.1 + (clearcoat_gloss)*0.001;
    Vector3 proj_vec = to_local(frame, half_vector);
    Real d_c = (pow(alpha_g,2) - 1)/(c_PI*log(pow(alpha_g,2))*(1 + (pow(alpha_g,2) - 1)*(pow(proj_vec.z,2))));

    Vector3 local_in = to_local(frame, dir_in);
    Vector3 local_out = to_local(frame, dir_out);

    Real delta_numerator_out = pow((1 + (pow(local_out.x*0.25,2) + pow(local_out.y*0.25,2))/(pow(local_out.z,2))),0.5) - 1.0;
    Real delta_numerator_in = pow((1 + (pow(local_in.x*0.25,2) + pow(local_in.y*0.25,2))/(pow(local_in.z,2))),0.5) - 1.0;
    Real delta_out = delta_numerator_out/2;
    Real delta_in  = delta_numerator_in/2;
    Real g_c = (1/(1 + delta_in))*(1/(1 + delta_out));

    Real const_denom = 4*(fabs(n_dot_in));

    return make_const_spectrum((f_c*d_c*g_c)/(const_denom));
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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

    Real const_denom = 4*(fabs(n_dot_out));

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, 
                               vertex.uv, 
                               vertex.uv_screen_size, 
                               texture_pool);

    Real alpha_g = (1-clearcoat_gloss)*0.1 + (clearcoat_gloss)*0.001;
    
    Vector3 proj_vec = to_local(frame, half_vector);
    Real d_c = (pow(alpha_g,2) - 1)/(c_PI*log(pow(alpha_g,2))*(1 + (pow(alpha_g,2) - 1)*(pow(proj_vec.z,2))));

    return (d_c*fabs(n_dot_h))/const_denom;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, 
                               vertex.uv, 
                               vertex.uv_screen_size, 
                               texture_pool);

    Real alpha_g = (1-clearcoat_gloss)*0.1 + (clearcoat_gloss)*0.001;
    
    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal = sample_clearcoat(local_dir_in, alpha_g, rnd_param_uv);
    Vector3 half_vector = to_world(frame, local_micro_normal);
        
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
     return BSDFSampleRecord{
            reflected,
            Real(0) /* eta */, alpha_g /* roughness */
        };

}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}