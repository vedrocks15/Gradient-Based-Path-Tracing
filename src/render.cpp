#include "render.h"
#include "intersection.h"
#include "material.h"
#include "parallel.h"
#include "path_tracing.h"
#include "vol_path_tracing.h"
#include "pcg.h"
#include "progress_reporter.h"
#include "scene.h"
// Added : / you will need fftw3 http://www.fftw.org/ to compile
#include <fftw3.h>
#include <iostream>

/// Render auxiliary buffers e.g., depth.
Image3 aux_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    parallel_for([&](const Vector2i &tile) {
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Ray ray = sample_primary(scene.camera, Vector2((x + Real(0.5)) / w, (y + Real(0.5)) / h));
                RayDifferential ray_diff = init_ray_differential(w, h);
                if (std::optional<PathVertex> vertex = intersect(scene, ray, ray_diff)) {
                    Real dist = distance(vertex->position, ray.org);
                    Vector3 color{0, 0, 0};
                    if (scene.options.integrator == Integrator::Depth) {
                        color = Vector3{dist, dist, dist};
                    } else if (scene.options.integrator == Integrator::ShadingNormal) {
                        // color = (vertex->shading_frame.n + Vector3{1, 1, 1}) / Real(2);
                        color = vertex->shading_frame.n;
                    } else if (scene.options.integrator == Integrator::MeanCurvature) {
                        Real kappa = vertex->mean_curvature;
                        color = Vector3{kappa, kappa, kappa};
                    } else if (scene.options.integrator == Integrator::RayDifferential) {
                        color = Vector3{ray_diff.radius, ray_diff.spread, Real(0)};
                    } else if (scene.options.integrator == Integrator::MipmapLevel) {
                        const Material &mat = scene.materials[vertex->material_id];
                        const TextureSpectrum &texture = get_texture(mat);
                        auto *t = std::get_if<ImageTexture<Spectrum>>(&texture);
                        if (t != nullptr) {
                            const Mipmap3 &mipmap = get_img3(scene.texture_pool, t->texture_id);
                            Vector2 uv{modulo(vertex->uv[0] * t->uscale, Real(1)),
                                       modulo(vertex->uv[1] * t->vscale, Real(1))};
                            // ray_diff.radius stores approximatedly dpdx,
                            // but we want dudx -- we get it through
                            // dpdx / dpdu
                            Real footprint = vertex->uv_screen_size;
                            Real scaled_footprint = max(get_width(mipmap), get_height(mipmap)) *
                                                    max(t->uscale, t->vscale) * footprint;
                            Real level = log2(max(scaled_footprint, Real(1e-8f)));
                            color = Vector3{level, level, level};
                        }
                    }
                    img(x, y) = color;
                } else {
                    img(x, y) = Vector3{0, 0, 0};
                }
            }
        }
    }, Vector2i(num_tiles_x, num_tiles_y));

    return img;
}

Image3 path_render(const Scene &scene) {

    // Defining camera constants....
    int w = scene.camera.width, h = scene.camera.height;

    // Variable that will hold the final image...
    Image3 img(w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    // tqdm like helper...
    ProgressReporter reporter(num_tiles_x * num_tiles_y);

    // parallelisation (shared variable in thread : tile ??)
    parallel_for([&](const Vector2i &tile) {
        // Use a different rng stream for each thread. (some tile based looping)
        pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);

        // Looping across rows & cols...
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {

                // Final contribution holding variable for the primal image
                Spectrum radiance = make_zero_spectrum();
                int spp = 256;//scene.options.samples_per_pixel;

                // Looping over samples....
                for (int s = 0; s < spp; s++) {
                    radiance += path_tracing(scene, x, y, rng);
                }
                img(x, y) = radiance / Real(spp);
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
    return img;
}

Image3 vol_path_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    auto f = vol_path_tracing;
    if (scene.options.vol_path_version == 1) {
        f = vol_path_tracing_1;
    } else if (scene.options.vol_path_version == 2) {
        f = vol_path_tracing_2;
    } else if (scene.options.vol_path_version == 3) {
        f = vol_path_tracing_3;
    } else if (scene.options.vol_path_version == 4) {
        f = vol_path_tracing_4;
    } else if (scene.options.vol_path_version == 5) {
        f = vol_path_tracing_5;
    } else if (scene.options.vol_path_version == 6) {
        f = vol_path_tracing;
    }

    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    parallel_for([&](const Vector2i &tile) {
        // Use a different rng stream for each thread.
        pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Spectrum radiance = make_zero_spectrum();
                int spp = scene.options.samples_per_pixel;
                for (int s = 0; s < spp; s++) {
                    Spectrum L = f(scene, x, y, rng);
                    if (isfinite(L)) {
                        // Hacky: exclude NaNs in the rendering.
                        radiance += L;
                    }
                }
                img(x, y) = radiance / Real(spp);
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
    return img;
}

// GB-PT implementation...
// screened Poisson solver from http://grail.cs.washington.edu/projects/screenedPoissonEq/
void fourierSolve(int width, int height, 
        const double* imgData, const double* imgGradX, 
        const double* imgGradY, double dataCost,
        double* imgOut) {
    
    int nodeCount = width * height;
    double* fftBuff = (double*) fftw_malloc(sizeof(*fftBuff) * nodeCount);
    //compute two 1D lookup tables for computing the DCT of a 2D Laplacian on the fly
    double* ftLapY = (double*) fftw_malloc(sizeof(*ftLapY) * height);
    double* ftLapX = (double*) fftw_malloc(sizeof(*ftLapX) * width);
    for(int x = 0; x < width; x++) {
        ftLapX[x] = 2.0 * cos(M_PI * x / (width - 1));
    }
    for(int y = 0; y < height; y++) {
        ftLapY[y] = -4.0 + (2.0 * cos(M_PI * y / (height - 1)));
    }
    //Create a DCT-I plan for, which is its own inverse.
    fftw_plan fftPlan; 
    fftPlan = fftw_plan_r2r_2d(height, width, 
            fftBuff, fftBuff, 
            FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE); //use FFTW_PATIENT when plan can be reused
    
            for(int iChannel = 0; iChannel < 3; iChannel++) {
        int nodeAddr        = 0;
        int pixelAddr       = iChannel;
        int rightPixelAddr  = 3 + iChannel;
        int topPixelAddr    = (width * 3) + iChannel;
        double dcSum = 0.0;

        // compute h_hat from u, gx, gy (see equation 48 in Bhat's paper), as well as the DC term of u's DCT.
        for(int y = 0; y < height; y++)
            for(int x = 0; x < width;  x++, 
                    nodeAddr++, pixelAddr += 3, rightPixelAddr += 3, topPixelAddr += 3) {
                // Compute DC term of u's DCT without computing the whole DCT.
                double dcMult = 1.0;
                if((x > 0) && (x < width  - 1))
                    dcMult *= 2.0;
                if((y > 0) && (y < height - 1))
                    dcMult *= 2.0;
                dcSum += dcMult * imgData[pixelAddr];

                fftBuff[nodeAddr] = dataCost * imgData[pixelAddr];      

                // Subtract g^x_x and g^y_y, with boundary factor of -2.0 to account for boundary reflections implicit in the DCT
                if((x > 0) && (x < width - 1))
                    fftBuff[nodeAddr] -= (imgGradX[rightPixelAddr] - imgGradX[pixelAddr]);
                else
                    fftBuff[nodeAddr] -= (-2.0 * imgGradX[pixelAddr]);

                if((y > 0) && (y < height - 1))
                    fftBuff[nodeAddr] -= (imgGradY[topPixelAddr] - imgGradY[pixelAddr]);
                else
                    fftBuff[nodeAddr] -= (-2.0 * imgGradY[pixelAddr]);
            }
        //transform h_hat to H_hat by taking the DCT of h_hat
        fftw_execute(fftPlan);

        //compute F_hat using H_hat (see equation 29 in Bhat's paper)
        nodeAddr = 0;
        for(int y = 0; y < height; y++)
            for(int x = 0; x < width;  x++, nodeAddr++) {
                float ftLapResponse = ftLapY[y] + ftLapX[x]; 
                fftBuff[nodeAddr] /= (dataCost - ftLapResponse);
            }
        /* Set the DC term of the solution to the value computed above (i.e., the DC term of imgData). 
         * set dcSum to the desired average when dataCost=0
         */
        fftBuff[0] = dcSum;

        //transform F_hat to f_hat by taking the inverse DCT of F_hat
        fftw_execute(fftPlan);   
        double fftDenom = 4.0 * (width - 1) * (height - 1);
        pixelAddr = iChannel;
        for(int iNode = 0; iNode < nodeCount; iNode++, pixelAddr += 3) {
            imgOut[pixelAddr] = fftBuff[iNode] / fftDenom;  
        }
    }

    fftw_free(fftBuff);
    fftw_free(ftLapX);
    fftw_free(ftLapY);
    fftw_destroy_plan(fftPlan);
}


Image3 gradient_path_render(const Scene &scene) {

    // 1.] Defining camera constants : 
    int w = scene.camera.width, h = scene.camera.height;

    // 2.] Holder Image variables for primal image & estimated gradient image
    // Since each pixel has 4 neighbours thus we have 4 extra variables.
    Image3 img(w, h);
    Image3 cx0(w, h);
    Image3 cy0(w, h);
    Image3 cx1(w, h);
    Image3 cy1(w, h);

    // 3.] Tiling added for parallelization
    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;
    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    
    // parallelisation, looping across tiles
    parallel_for([&](const Vector2i &tile) {

        // Looping inside each tile.
        // Use a different rng stream for each thread.
        pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        
        // Tile start & end
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w); // important if tile overshoots the image dims
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);

        //4.] Looping across rows & cols.
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {

                int spp = 1000 ;//scene.options.samples_per_pixel;

                // Variables to summate the contribution from each sample per pixel.
                Spectrum r = make_zero_spectrum();
                Spectrum rdX0 = make_zero_spectrum();
                Spectrum rdX1 = make_zero_spectrum();
                Spectrum rdY0 = make_zero_spectrum();
                Spectrum rdY1 = make_zero_spectrum();

                // Looping over samples....
                for (int s = 0; s < spp; s++) {

                    // Logic for original path tracing & shift mapping...
                    // We are computing all inside one function instead 
                    // of first generating the base path, followed by computing the shift path
                    // & then the jacobians + gradients.
                    GraidentPTRadiance holder_var = grad_path_tracing(scene, x, y, rng);

                    if (holder_var.prob > 0.0)
                    {
                        r = r + holder_var.radiance/Real(spp);
                        rdX0 = rdX0 + (holder_var.contrib - holder_var.contribX0) * (holder_var.wX0 / (holder_var.prob * Real(spp)));
                        rdY0 = rdY0 + (holder_var.contrib - holder_var.contribY0) * (holder_var.wY0 / (holder_var.prob * Real(spp)));
                        rdX1 = rdX1 + (holder_var.contribX1 - holder_var.contrib) * (holder_var.wX1 / (holder_var.prob * Real(spp)));
                        rdY1 = rdY1 + (holder_var.contribY1 - holder_var.contrib) * (holder_var.wY1 / (holder_var.prob * Real(spp)));
                    }
                
                }
                // Collecting all image & gradient combinations...
                // If path was not generated only a zero spectrum vector is added....
                img(x, y) += r;
                cx0(x, y) += rdX0;
                cy0(x, y) += rdY0;
                cx1(x, y) += rdX1;
                cy1(x, y) += rdY1;
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();


    // add possion logic ......
    Vector3 *cx = new Vector3[w*h], *cy = new Vector3[w*h];
    Vector3 *c = new Vector3[w*h];
    Vector3 *out = new Vector3[w*h];

    for (int y=0; y<h; y++)
        for (int x=0; x<w; x++) {
            int i = y * w + x;
            // flattening image
            c[i] = img(x,y);
            if (x == 0) cx[i] = cx0(x,y);
            else cx[i] = cx0(x,y) + cx1(x-1,y);
            
            if ( y== 0)cy[i] = cy0(x,y);
            else cy[i] = cy0(x,y) + cy1(x,y-1);
        }
    

    fourierSolve(w, h, (double*)c, (double*)cx, (double*)cy, 0.04, (double*)out);

    Image3 fout(w, h); 
    Image3 fc(w, h); 
    Image3 fcx(w, h); 
    Image3 fcy(w, h); 

    for (int y=0; y<h; y++)
        for (int x=0; x<w; x++) {
            int i = y * w + x;
            fout(x,y) = out[i];
            fc(x,y) = c[i];
            fcx(x,y).x = fabs(cx[i].x); fcx(x,y).y = fabs(cx[i].y); fcx(x,y).z = fabs(cx[i].z);
            fcy(x,y).x = fabs(cy[i].x); fcy(x,y).y = fabs(cy[i].y); fcy(x,y).z = fabs(cy[i].z);
        }

    return fout;
}



Image3 render(const Scene &scene) {
    if (scene.options.integrator == Integrator::Depth ||
            scene.options.integrator == Integrator::ShadingNormal ||
            scene.options.integrator == Integrator::MeanCurvature ||
            scene.options.integrator == Integrator::RayDifferential ||
            scene.options.integrator == Integrator::MipmapLevel) {
        return aux_render(scene);
    } else if (scene.options.integrator == Integrator::Path) {
        return path_render(scene); 
    } else if (scene.options.integrator == Integrator::VolPath) {
        return vol_path_render(scene);
    }
    // // New Implementation added for gradient based rendering....
    else if(scene.options.integrator == Integrator::GradPath){
        return gradient_path_render(scene);
    } 
    else {
        assert(false);
        return Image3();
    }
}
