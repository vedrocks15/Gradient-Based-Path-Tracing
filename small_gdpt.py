
// smallgdpt: a simple implementation of gradient domain path tracing 
//                                       https://mediatech.aalto.fi/publications/graphics/GPT/ 
// adapted from smallpt by Kevin Beason http://www.kevinbeason.com/smallpt/
// and a screened poisson solver by Pravin Bhat http://grail.cs.washington.edu/projects/screenedPoissonEq/
// to build, type: g++ -o smallgdpt -fopenmp -O3 smallgdpt.cpp -L/usr/local/lib -lm -lfftw3
// you will need fftw3 http://www.fftw.org/ to compile
// usage: ./smallgdpt [number of samples per pixel]
#include <fftw3.h>
#include <math.h>   
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>  

const int MAX_DEPTH = 32; 

struct Vec {        
    double x, y, z;                  // position, also color (r,g,b)
    Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
    Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
    Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
    Vec operator-() const {return Vec(-x,-y,-z);}
    Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
    Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
    Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
    double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross:
    Vec operator%(Vec&b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
    double max() const {return x>y && x>z ? x : y > z ? y : z;}
};

struct Ray { Vec o, d; Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };

enum Refl_t { DIFF, SPEC, REFR };  // material types

struct Sphere {
    double rad;       // radius
    Vec p, e, c;      // position, emission, color
    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
        rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &r) const { // returns distance, 0 if nohit
        Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
        if (det<0) return 0; else det=sqrt(det);
        return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
    }
};

int width = 1024; int height = 768;
Ray cam(Vec(50,50,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
Vec cx=Vec(width*.5135/height), cy=(cx%cam.d).norm()*.5135;
Sphere spheres[] = {//Scene: radius, position, emission, color, material
    Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
    Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
    Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
    Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
    Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
    Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1.0,1.0,1.0)*.999, SPEC),//Mirr
    Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1.0,1.0,1.0)*.999, REFR),//Glas
    Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite
};

struct PathVert {
    Vec p; Vec n; int id;
};

struct Path {
    PathVert verts[MAX_DEPTH];
    double rnds[2*MAX_DEPTH];
    int vertCount;
    int x, y;
};

inline bool intersect(const Ray &r, double &t, int &id){
    double n=sizeof(spheres)/sizeof(Sphere), d, inf=t=1e20;
    for(int i=int(n);i--;) if((d=spheres[i].intersect(r))&&d<t){t=d;id=i;}
    return t<inf;
}

Vec reflect(const Vec &d, const Vec &n) {
    return d - n * 2.0 * n.dot(d);
}

Ray sampleBSDF(const Ray &ray, const Sphere &obj, const PathVert &vert, double u0, double u1) {
    if (obj.refl == DIFF) {  
        double r1=2*M_PI*u0, r2=u1, r2s=sqrt(r2);
        Vec nl=vert.n.dot(ray.d)<0?vert.n:vert.n*-1; // flip normal if needed
        Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
        Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
        return Ray(vert.p, d);      
    } else if (obj.refl == SPEC) {
        return Ray(vert.p, reflect(ray.d, vert.n));
    } else { //REFR
        Ray reflRay(vert.p, reflect(ray.d, vert.n));              
        bool into = vert.n.dot(ray.d)<0;
        Vec nl = into?vert.n:vert.n*-1;
        double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=ray.d.dot(nl), cos2t;
        if ((cos2t = 1-nnt*nnt*(1-ddn*ddn)) < 0) { // total internal reflection
            return reflRay;
        }
        Vec tdir = (ray.d*nnt - vert.n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
        double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(vert.n));
        double Re=R0+(1-R0)*c*c*c*c*c, P=.25+.5*Re; // schlick
        if (u0 < P) {
            return reflRay;
        } else {
            return Ray(vert.p, tdir);
        }
    }
    return Ray(Vec(), Vec());
}

double BSDFProb(const Refl_t &refl, const Vec &wi, const Vec &n, const Vec &wo) {
    if (refl == DIFF) {
        double cosTheta = fabs(wo.dot(n));
        return (cosTheta/M_PI);
    } else if (refl == SPEC) {
        return 1.0;
    } else { //REFR
        bool refl = wi.dot(n) * wo.dot(n) > 0.0;
        bool into = n.dot(wi) > 0;
        Vec nl = into ? n : n*-1; // flip normal if needed
        Vec d = -wi;
        double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=d.dot(nl), cos2t;
        double P = refl ? 1.0 : 0.0;
        if ((cos2t = 1-nnt*nnt*(1-ddn*ddn)) > 0) {
            Vec tdir = (d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
            double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
            double Re=R0+(1-R0)*c*c*c*c*c;
            P = .25+.5*Re; if (!refl) P = 1.0 - P;            
        }
        return P;
    }
    return 0.0;
}

// generate a light path from scratch
bool generatePath(int x, int y, unsigned short *rng, Path &path) {  
    path.x = x; path.y = y;
    path.rnds[0] = erand48(rng); path.rnds[1] = erand48(rng); 
    Vec d = cx*( (path.rnds[0] + x)/width - .5) +
            cy*( (path.rnds[1] + y)/height - .5) + cam.d;  
    // Camera rays are pushed forward to start in interior  
    Ray ray(cam.o+d*140, d.norm());
    path.vertCount = 0;
    for (int depth = 1; depth <= MAX_DEPTH; depth++) {
        double t; int id = -1;
        if (!intersect(ray, t, id)) return false;    
        const Sphere &obj = spheres[id];        
        PathVert vert; 
        vert.p = ray.o + ray.d*t; vert.n = (vert.p - obj.p).norm(); vert.id = id;    
        path.verts[depth-1] = vert; path.vertCount++;    
        double p = obj.c.max();    
        if (p <= 0.0) return true; // assume refl=0 -> light source
        if (depth == MAX_DEPTH) return false;
        path.rnds[2*depth] = erand48(rng); path.rnds[2*depth+1] = erand48(rng);
        ray = sampleBSDF(ray, obj, vert, path.rnds[2*depth], path.rnds[2*depth+1]);
    }
    return false;
}

// "shift" a light path to a specific pixel
bool shiftPath(int x, int y, const Path &basePath, Path &offsetPath, double &jacobian) {    
    offsetPath.x = x; offsetPath.y = y;
    Vec baseWi = -(cx*( (basePath.rnds[0] + basePath.x)/width - .5) +
                   cy*( (basePath.rnds[1] + basePath.y)/height - .5) + cam.d).norm();    
    Vec d = cx*( (basePath.rnds[0] + x)/width - .5) +
            cy*( (basePath.rnds[1] + y)/height - .5) + cam.d;    
    Ray ray(cam.o+d*140, d.norm());
    Vec wi = -ray.d;
    offsetPath.vertCount = basePath.vertCount;  
    memcpy(offsetPath.verts, basePath.verts, sizeof(PathVert) * basePath.vertCount);
    jacobian = 1.0;
    for (int vertId = 0; vertId < basePath.vertCount; vertId++) {
        int depth = vertId + 1;
        double t; int id = -1;
        if (!intersect(ray, t, id)) return false;        
        const Sphere &obj = spheres[id]; 
        const Sphere &baseObj = spheres[basePath.verts[vertId].id];
        if (obj.refl != baseObj.refl) return false;
        PathVert vert; 
        vert.p = ray.o + ray.d*t; vert.n = (vert.p - obj.p).norm(); vert.id = id;    
        offsetPath.verts[vertId] = vert; 
        if (vertId == basePath.vertCount - 1) break;
        if (obj.refl == DIFF && spheres[basePath.verts[vertId + 1].id].refl == DIFF) {
            // connect back to base path, jacobian = ratio of geometry term
            if (!intersect(Ray(vert.p, (basePath.verts[depth].p - vert.p).norm()), t, id) || 
                    id != basePath.verts[vertId + 1].id) return false;        

            Vec baseP0 = basePath.verts[depth - 1].p;
            Vec p1 = basePath.verts[depth].p;
            Vec baseN0 = basePath.verts[depth - 1].n;
            Vec n1 = basePath.verts[depth].n;
            Vec baseDir = p1 - baseP0;
            double baseDist2 = baseDir.dot(baseDir);
            baseDir = baseDir * (1.0 / sqrt(baseDist2));
            double baseGeom = fabs(baseDir.dot(n1)) * fabs(baseDir.dot(baseN0)) / baseDist2;
            Vec shiftDir = p1 - vert.p;
            double shiftDist2 = shiftDir.dot(shiftDir);     
            shiftDir = shiftDir * (1.0 / sqrt(shiftDist2));
            double shiftGeom = fabs(shiftDir.dot(n1)) * fabs(shiftDir.dot(vert.n)) / shiftDist2;
            jacobian *= (shiftGeom / baseGeom);
            return true;
        }
        
        // copy the random numbers used to sample BRDF, jacobian = ratio of inverse PDF
        // this should be simpler than the half-vector based shift described in the paper
        ray = sampleBSDF(ray, obj, vert, basePath.rnds[2*depth], basePath.rnds[2*depth+1]);
        Vec baseWo = (basePath.verts[vertId + 1].p - basePath.verts[vertId].p).norm();
        double basePDF = BSDFProb(baseObj.refl, baseWi, basePath.verts[vertId].n, baseWo);
        double shiftPDF = BSDFProb(obj.refl, wi, vert.n, ray.d);
        if (shiftPDF <= 0.0) return false;
        jacobian *= (basePDF / shiftPDF);
        baseWi = -baseWo; wi = -ray.d;
    }   
    const Sphere &obj = spheres[offsetPath.verts[offsetPath.vertCount-1].id];
    double p = obj.c.max();    
    return p <= 0.0; // assume refl=0 -> light source
}

// path contribution in solid angle domain
Vec pathContrib(const Path &path) {   
    Vec throughput(1,1,1);
    Vec wi = -(cx*( (path.rnds[0] + path.x)/width - .5) +
               cy*( (path.rnds[1] + path.y)/height - .5) + cam.d).norm();    
    for (int vert = 0; vert < path.vertCount - 1; vert++) {
        const PathVert &currVert = path.verts[vert];
        const PathVert &nextVert = path.verts[vert + 1];
        Vec wo = (nextVert.p - currVert.p).norm();
        double cosTheta = fabs(wo.dot(currVert.n));
        const Sphere &obj = spheres[path.verts[vert].id];     
        if (cosTheta <= 1e-6) return Vec();
        if (obj.refl == DIFF) {
            throughput = throughput.mult(obj.c*(cosTheta/M_PI)); 
        } else if (obj.refl == SPEC) {
            throughput = throughput.mult(obj.c);
        } else { //REFR            
            bool refl = wi.dot(currVert.n) * wo.dot(currVert.n) > 0.0;
            bool into = currVert.n.dot(wi) > 0;
            Vec d = -wi;
            Vec nl = into ? currVert.n : currVert.n*-1; // flip normal if needed
            double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=d.dot(nl), cos2t;
            double fresnel = refl ? 1.0 : 0.0;
            if ((cos2t = 1-nnt*nnt*(1-ddn*ddn)) > 0) {
                Vec tdir = (d*nnt - currVert.n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
                double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(currVert.n));
                double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re;              
                fresnel = refl ? Re : Tr;
            }
            throughput = throughput.mult(obj.c * fresnel); 
        }
        wi = -wo;
    }    
    const Sphere &obj = spheres[path.verts[path.vertCount-1].id];
    return throughput.mult(obj.e);  
}

// path probability in solid angle domain
double pathProb(const Path &path) {
    Vec wi = -(cx*( (path.rnds[0] + path.x)/width - .5) +
               cy*( (path.rnds[1] + path.y)/height - .5) + cam.d).norm();    
    double prob = 1.0;
    for (int vert = 0; vert < path.vertCount - 1; vert++) {    
        const PathVert &currVert = path.verts[vert];
        const PathVert &nextVert = path.verts[vert + 1];
        Vec wo = (nextVert.p - currVert.p).norm();
        double cosTheta = fabs(wo.dot(currVert.n));
        const Sphere &obj = spheres[path.verts[vert].id];     
        if (cosTheta <= 1e-6) return 0.0;       
        prob *= BSDFProb(obj.refl, wi, currVert.n, wo);
        if (prob <= 0.0) return 0.0;
        wi = -wo;
    }
    return prob;
}

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

int main(int argc, char *argv[]){
    int samps = argc==2 ? atoi(argv[1]) : 4; // # samples
    Vec *c=new Vec[width * height];
    Vec *cx0=new Vec[width * height];
    Vec *cy0=new Vec[width * height];
    Vec *cx1=new Vec[width * height];
    Vec *cy1=new Vec[width * height];
#pragma omp parallel for schedule(dynamic, 1) // OpenMP
    for (int y=0; y<height; y++){                       // Loop over image rows
        fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps,100.*y/(height-1));
        for (unsigned short x=0, rng[3]={0,0,y*y*y}; x<width; x++) {  // Loop cols
            Vec r, rdx0, rdy0, rdx1, rdy1;
            for (int s=0; s<samps; s++){
                Path path, oPath; double jacobian;        
                if(generatePath(x, y, rng, path)) {                     
                    Vec contrib = pathContrib(path);
                    double prob = pathProb(path);
                    if (prob > 0.0) {                              
                        Vec contribX0, contribY0;
                        Vec contribX1, contribY1;
                        double wX0 = 1, wY0 = 1;
                        double wX1 = 1, wY1 = 1; 
                        r = r + (contrib * (1.0 / prob)) * (1.0 / (double)samps);          
                        if(shiftPath(x-1, y, path, oPath, jacobian)) {
                            contribX0 = pathContrib(oPath) * jacobian;
                            double pX0 = pathProb(oPath) * jacobian;
                            wX0 = prob / (prob + pX0);
                        }
                        if(shiftPath(x, y+1, path, oPath, jacobian)) {
                            contribY0 = pathContrib(oPath) * jacobian;
                            double pY0 = pathProb(oPath) * jacobian;
                            wY0 = prob / (prob + pY0);
                        }
                        if(shiftPath(x+1, y, path, oPath, jacobian)) {
                            contribX1 = pathContrib(oPath) * jacobian;
                            double pX1 = pathProb(oPath) * jacobian;
                            wX1 = prob / (prob + pX1);
                        }
                        if(shiftPath(x, y-1, path, oPath, jacobian)) {
                            contribY1 = pathContrib(oPath) * jacobian;
                            double pY1 = pathProb(oPath) * jacobian;
                            wY1 = prob / (prob + pY1);
                        }                       

                        rdx0 = rdx0 + (contrib - contribX0) * (wX0 / (prob * (double)samps));
                        rdy0 = rdy0 + (contrib - contribY0) * (wY0 / (prob * (double)samps));
                        rdx1 = rdx1 + (contribX1 - contrib) * (wX1 / (prob * (double)samps));
                        rdy1 = rdy1 + (contribY1 - contrib) * (wY1 / (prob * (double)samps));
                    }
                }
            }                   
            int i = (height - y - 1) * width + x;
            c[i]  = c[i] + r; 
            cx0[i] = cx0[i] + rdx0;  cy0[i] = cy0[i] + rdy0;  
            cx1[i] = cx1[i] + rdx1;  cy1[i] = cy1[i] + rdy1;
        }
    }
    Vec *cx=new Vec[width * height], *cy = new Vec[width * height];
    for (int y=0; y<height; y++)
        for (int x=0; x<width; x++) {
            int i = y * width + x;
            if (x == 0) cx[i] = cx0[i];
            else cx[i] = cx0[i] + cx1[i-1];
            if (y == 0) cy[i] = cy0[i];
            else cy[i] = cy0[i] + cy1[i-width];
        }
    Vec *out=new Vec[width * height];
    fourierSolve(width, height, (double*)c, (double*)cx, (double*)cy, 0.04, (double*)out);

    int npixel = 3 * width * height;
    float *fc   = new float[npixel], *fout = new float[npixel];
    float *fcx  = new float[npixel], *fcy  = new float[npixel];
    for(int i = 0; i < width * height; i++) { //pfm requires single precision
        fc[3*i]   = c[i].x;   fc[3*i+1]   = c[i].y;   fc[3*i+2]   = c[i].z;    
        fout[3*i] = out[i].x; fout[3*i+1] = out[i].y; fout[3*i+2] = out[i].z;
        fcx[3*i] = fabs(cx[i].x); fcx[3*i+1] = fabs(cx[i].y); fcx[3*i+2] = fabs(cx[i].z);
        fcy[3*i] = fabs(cy[i].x); fcy[3*i+1] = fabs(cy[i].y); fcy[3*i+2] = fabs(cy[i].z);
    }
    FILE *f = fopen("image.pfm", "w");         // Write image to PFM files.
    fprintf(f, "PF\n%d %d\n%d\n", width, height, -1);
    fwrite(fc, sizeof(float), npixel, f); fclose(f);
    f = fopen("image_dx.pfm", "w");         
    fprintf(f, "PF\n%d %d\n%d\n", width, height, -1);
    fwrite(fcx, sizeof(float), npixel, f); fclose(f);
    f = fopen("image_dy.pfm", "w");         
    fprintf(f, "PF\n%d %d\n%d\n", width, height, -1);
    fwrite(fcy, sizeof(float), npixel, f); fclose(f); 
    f = fopen("image_poisson.pfm", "w");
    fprintf(f, "PF\n%d %d\n%d\n", width, height, -1);
    fwrite(fout, sizeof(float), npixel, f); fclose(f);  
    return 0;
}