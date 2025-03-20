#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1....
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Spectrum base_color_comp = eval(bsdf.base_color, 
                               vertex.uv, 
                               vertex.uv_screen_size, 
                               texture_pool);

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real anisotropic = eval(bsdf.anisotropic, 
                          vertex.uv, 
                          vertex.uv_screen_size, 
                          texture_pool);

    // Compute F / D / G
    // Note that we use the incoming direction
    // for evaluating the Fresnel reflection amount.
    // We can also use outgoing direction -- then we would need to
    // use 1/bsdf.eta and we will get the same result.
    // However, using the incoming direction allows
    // us to use F to decide whether to reflect or refract during sampling.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    // X_proj :
    Vector3 proj_vec = to_local(frame, half_vector);
    
    Real aspect  = pow((1 - 0.9*anisotropic),0.5);
    Real alpha_x = max(0.0001, pow(roughness,2)/aspect); 
    Real alpha_y = max(0.0001, pow(roughness,2)*aspect); 

    Real dist_const = c_PI*alpha_x*alpha_y;
    Real d_m = 1 /(dist_const*pow((pow(proj_vec.x,2)/(pow(alpha_x,2)) +  pow(proj_vec.y,2)/(pow(alpha_y,2)) + pow(proj_vec.z,2)),2));

    Vector3 local_in = to_local(frame, dir_in);
    Vector3 local_out = to_local(frame, dir_out);

    Real internal_out = (pow(local_out.x*alpha_x,2) + pow(local_out.y*alpha_y,2)) / (pow(local_out.z,2));
    Real internal_in  = (pow(local_in.x*alpha_x,2) + pow(local_in.y*alpha_y,2)) / (pow(local_in.z,2));

    Real delta_out = (pow(1 + internal_out,0.5)-1)/2;
    Real delta_in  = (pow(1 + internal_in,0.5)-1)/2;

    Real g_m = (1/(1 + delta_in))*(1/(1 + delta_out));


    
    if (reflect) {
        return base_color_comp * (F * d_m * g_m) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real h_dot_out = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Spectrum color_sq = make_zero_spectrum();
        color_sq.x = pow(base_color_comp.x,0.5);
        color_sq.y = pow(base_color_comp.y,0.5);
        color_sq.z = pow(base_color_comp.z,0.5);
        return  color_sq*((1 - F) * d_m * g_m * fabs(h_dot_out * h_dot_in)) / 
            (fabs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom);
    }


}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    
    Vector3 proj_vec = to_local(frame, half_vector);

     Real anisotropic = eval(bsdf.anisotropic, 
                          vertex.uv, 
                          vertex.uv_screen_size, 
                          texture_pool);
    
    Real aspect  = pow((1 - 0.9*anisotropic),0.5);
    Real alpha_x = max(0.0001, pow(roughness,2)/aspect); 
    Real alpha_y = max(0.0001, pow(roughness,2)*aspect); 

    Real dist_const = c_PI*alpha_x*alpha_y;
    Real d_m = 1 /(dist_const*pow((pow(proj_vec.x,2)/(pow(alpha_x,2)) +  pow(proj_vec.y,2)/(pow(alpha_y,2)) + pow(proj_vec.z,2)),2));

    Vector3 local_in = to_local(frame, dir_in);
    Vector3 local_out = to_local(frame, dir_out);

    Real internal_out = (pow(local_out.x*alpha_x,2) + pow(local_out.y*alpha_y,2)) / (pow(local_out.z,2));
    Real internal_in  = (pow(local_in.x*alpha_x,2) + pow(local_in.y*alpha_y,2)) / (pow(local_in.z,2));

    Real delta_out = (pow(1 + internal_out,0.5)-1)/2;
    Real delta_in  = (pow(1 + internal_in,0.5)-1)/2;

    Real g_m = (1/(1 + delta_in));



    if (reflect) {
        return (F * d_m * g_m) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real h_dot_out = dot(half_vector, dir_out);
        Real h_dot_in = dot(half_vector, dir_in);
        Real prod_val = fabs(h_dot_out*h_dot_in);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real frame_in = dot(frame.n, dir_in);

        //Real sqrt_denom = h_dot_in + eta * h_dot_out;
        //Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return ((1 - F) * d_m * g_m * fabs(prod_val))/(fabs(frame_in) * sqrt_denom*sqrt_denom);
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Real alpha = roughness * roughness;
    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    if (rnd_param_w <= F) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
    } else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            return {};
        }
        // flip half_vector if needed
        if (h_dot_in < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_out= sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        return BSDFSampleRecord{refracted, eta, roughness};
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
