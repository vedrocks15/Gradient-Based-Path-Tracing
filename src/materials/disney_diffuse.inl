Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    
   
    Real roughness = eval(bsdf.roughness, 
                          vertex.uv, 
                          vertex.uv_screen_size, 
                          texture_pool);
    
    Spectrum base_color_comp = eval(bsdf.base_color, 
                               vertex.uv, 
                               vertex.uv_screen_size, 
                               texture_pool);


    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // Base diffuse component.....
    Real f_d_90  = 0.5 + 2*roughness*pow(fabs(half_out),2);
    Real const_pow_5_out = pow((1 -  fabs(n_dot_out)),5);
    Real const_pow_5_in  = pow((1 -  fabs(n_dot_in)),5);

    Real f_d_out = 1.0 + (f_d_90 - 1.0)*const_pow_5_out;
    Real f_d_in  = 1.0 + (f_d_90 - 1.0)*const_pow_5_in;

    Spectrum f_base_diffuse = (base_color_comp*f_d_in*f_d_out*fabs(n_dot_out))/c_PI;

    // Sub-surface component.....
    Real f_ss_90 = roughness*pow(fabs(half_out),2);
    Real f_ss_in  = 1.0 + (f_ss_90 - 1.0)*const_pow_5_in;
    Real f_ss_out = 1.0 + (f_ss_90 - 1.0)*const_pow_5_out;

    Spectrum f_sub_surface = (1.25*(base_color_comp)/c_PI)*((f_ss_in*f_ss_out)*((1/(fabs(n_dot_in) + fabs(n_dot_out))) - make_const_spectrum(0.5)) +  make_const_spectrum(0.5))*fabs(n_dot_out);
    
    Real subsurface_val = eval(bsdf.subsurface, 
                            vertex.uv, 
                            vertex.uv_screen_size, 
                            texture_pool);


    return ((1-subsurface_val)*f_base_diffuse + subsurface_val*f_sub_surface);
    
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
    
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Real roughness = eval(bsdf.roughness, 
                          vertex.uv, 
                          vertex.uv_screen_size, 
                          texture_pool);

    
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // // Homework 1: implement this! (cosine weighted)
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, roughness /* roughness */};
    
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
