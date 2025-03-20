#include "../microfacet.h"


Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    bool reflect_new_condition = dot(vertex.geometric_normal, dir_in) <= 0;

    //Get individual components...
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic  = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular  = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum c_tint = make_const_spectrum(1.0);
    if (luminance(base_color) > 0)
    {
        c_tint = base_color/luminance(base_color);
    }

    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Spectrum K_s = (1- specular_tint) + specular_tint * c_tint;

    Real r_0 = pow((eta-1), 2)/pow((eta+1),2);
    Spectrum c_0 = specular * r_0 * (1-metallic)* K_s + metallic * base_color;

    // Get disney components :
    Spectrum disney_diffuse_spectrum = make_const_spectrum(0);
    Spectrum disney_metallic_spectrum = make_const_spectrum(0);
    DisneyGlass disney_glass = DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, eta};
    Spectrum disney_glass_spectrum = eval_op::operator()(disney_glass);
    Spectrum disney_clearcoat_spectrum = make_const_spectrum(0);
    Spectrum disney_sheen_spectrum = make_const_spectrum(0);

    Real diffuseWeight   = (1 - specular_transmission)*(1 - metallic);
    Real metalWeight     = (1 - specular_transmission*(1 - metallic));
    Real clearcoatWeight =  0.25*clearcoat;  
    Real glassWeight     = (1 - metallic)*specular_transmission;
    Real sheenWeight     = (1-metallic)*sheen;

    if (reflect_new_condition)
    {
        return glassWeight*disney_glass_spectrum;
    }

    DisneyDiffuse disney_diffuse = DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    DisneyMetal disney_metallic = DisneyMetal{make_constant_spectrum_texture(c_0), bsdf.roughness, bsdf.anisotropic};
    DisneyClearcoat disney_clearcoat = DisneyClearcoat{bsdf.clearcoat_gloss};
    DisneySheen disney_sheen = DisneySheen{bsdf.base_color, bsdf.sheen_tint};
       
    disney_diffuse_spectrum = eval_op::operator()(disney_diffuse);
    disney_metallic_spectrum = eval_op::operator()(disney_metallic);
    disney_sheen_spectrum = eval_op::operator()(disney_sheen);
    disney_clearcoat_spectrum = eval_op::operator()(disney_clearcoat);

   
    return diffuseWeight*disney_diffuse_spectrum + 
        metalWeight*disney_metallic_spectrum+
        clearcoatWeight*disney_clearcoat_spectrum+
        glassWeight*disney_glass_spectrum+
        sheenWeight*disney_sheen_spectrum;

}


Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1: implement this!
    bool reflect_new_condition = dot(vertex.geometric_normal, dir_in) <= 0;

    //Get individual components...
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic  = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular  = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum c_tint = make_const_spectrum(1.0);
    if (luminance(base_color) > 0)
    {
        c_tint = base_color/luminance(base_color);
    }

    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Spectrum K_s = (1- specular_tint) + specular_tint * c_tint;

    Real r_0 = pow((eta-1), 2)/pow((eta+1),2);
    Spectrum c_0 = specular * r_0 * (1-metallic)* K_s + metallic * base_color;

    // Get disney components :
    Real disney_diffuse_real  = 0.0;
    Real disney_metallic_real = 0.0;
    Real disney_glass_real = 0.0;
    Real disney_clearcoat_glass = 0.0;

    Real diffuseWeight   = (1 - specular_transmission)*(1 - metallic);
    Real metalWeight     = (1 - specular_transmission*(1 - metallic));
    Real clearcoatWeight =  0.25*clearcoat;  
    Real glassWeight     = (1 - metallic)*specular_transmission;

    Real net_weight = diffuseWeight + metalWeight + clearcoatWeight+ glassWeight;

    DisneyGlass disney_glass = DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, eta};

    if (reflect_new_condition)
    {
        // Normalised....
        return pdf_sample_bsdf_op::operator()(disney_glass); 
    }

    DisneyDiffuse disney_diffuse = DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    DisneyMetal disney_metallic = DisneyMetal{make_constant_spectrum_texture(c_0), bsdf.roughness, bsdf.anisotropic};
    DisneyClearcoat disney_clearcoat = DisneyClearcoat{bsdf.clearcoat_gloss};
       
   
    return (diffuseWeight/net_weight)*pdf_sample_bsdf_op::operator()(disney_diffuse) +
           (metalWeight/net_weight)*pdf_sample_bsdf_op::operator()(disney_metallic) +
           (clearcoatWeight / net_weight)*pdf_sample_bsdf_op::operator()(disney_clearcoat);
           (glassWeight / net_weight)*pdf_sample_bsdf_op::operator()(disney_glass);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }


    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic  = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular  = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum c_tint = make_const_spectrum(1.0);
    if (luminance(base_color) > 0)
    {
        c_tint = base_color/luminance(base_color);
    }


    Real eta =  dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Spectrum K_s = (1- specular_tint) + specular_tint * c_tint;
    Real r_0 = pow((eta-1), 2)/pow((eta+1),2);
    Spectrum c_0 = specular * r_0 * (1-metallic)* K_s + metallic * base_color;

    DisneyDiffuse disney_diffuse = DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    DisneyMetal disney_metallic = DisneyMetal{make_constant_spectrum_texture(c_0), bsdf.roughness, bsdf.anisotropic};
    DisneyClearcoat disney_clearcoat = DisneyClearcoat{bsdf.clearcoat_gloss};
    DisneyGlass disney_glass = DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, eta};

    // // Homework 1: implement this!
    // pcg32_state rng_1 = init_pcg32();
    Real rand_num = rnd_param_uv[0];
    
    // Diffuse : 
    if (rand_num < 0.25)
    {
        return sample_bsdf_op::operator()(disney_diffuse);
    }
    else
    if((rand_num >= 0.25) && (rand_num < 0.5))
    {
        return sample_bsdf_op::operator()(disney_metallic);
    }
    else
    if((rand_num >= 0.5) && (rand_num < 0.75))
    {
        return sample_bsdf_op::operator()(disney_clearcoat);
    }
    //return sample_bsdf_op::operator()(disney_glass);
    return sample_bsdf_op::operator()(disney_clearcoat);

}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
