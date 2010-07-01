// includes from OpenEngine base
#include <Logging/Logger.h>
#include <Logging/StreamLogger.h>
#include <Resources/Texture2D.h>
#include <Resources/Tex.h>
#include <Resources/ResourceManager.h>
#include <Resources/DirectoryManager.h>
#include <Utils/Timer.h>
#include <Utils/TexUtils.h>

// includes from OpenEngine extensions
#include <Resources/FreeImage.h>
#include <Utils/TextureTool.h>

// other includes
#include <stdlib.h>

using namespace OpenEngine;
using namespace Logging;

//typedef double REAL;
typedef float REAL;

REAL abso(REAL n) {
    if (n < 0) n *= -1; 
    return n;
}

REAL HeatEquation(REAL dt, FloatTexture2DPtr ut) {
    unsigned int w = ut->GetWidth();
    unsigned int h = ut->GetHeight();
    FloatTexture2DPtr u0(ut->Clone());
    
    REAL error = 0.0f;
    for (unsigned int x=1; x<w-1; x++) {
        for (unsigned int y=1; y<h-1; y++) {
            REAL u0_xx = (*u0)(x+1,y) - 2*(*u0)(x,y) + (*u0)(x-1,y);
            REAL u0_yy = (*u0)(x,y+1) - 2*(*u0)(x,y) + (*u0)(x,y-1);
            REAL Lu0 = u0_xx + u0_yy;
            error += abso(Lu0);
            // #warning HEATEQ: CFG condition on dt!
            (*ut)(x,y) = (*u0)(x,y) + Lu0*dt;
        }
    }
    
    // solve boundary conditions
    for (unsigned int y=1; y<h-1; y++) {
        (*ut)(0,y) = (*ut)(1,y);
        (*ut)(w-1,y) = (*ut)(w-2,y);
    }
    for (unsigned int x=0; x<w; x++) {
        (*ut)(x,0) = (*ut)(x,1);
        (*ut)(x,h-1) = (*ut)(x,h-2);
    }

    //delete u0;
    return error;
}

FloatTexture2DPtr J(REAL dt, FloatTexture2DPtr input) {
    FloatTexture2DPtr calc(input->Clone());

    static unsigned int count = 0;
    static bool done = false;
    static REAL lastError = 0;
    //for (unsigned int i=0; i<1000;i++) {
    //std::ofstream file;
    //file.open("graph.gnu");
    
    while (!done) {
        REAL error = HeatEquation(dt, calc);
        REAL dError = error-lastError;
        //file << count << " " << error << std::endl;
        logger.info << "Heat Equation (iteration " << count++ << "), ";
        logger.info << "error: " << error << ", ";
        logger.info << "delta error: " << dError << logger.end;
        if (abso(error) < 250) done = true;
        lastError = error;
    }
    //file.close();
    
    return calc;
}

FloatTexture2DPtr ROFIteration(REAL dt,
               FloatTexture2DPtr u0,
               FloatTexture2DPtr v0,
               FloatTexture2DPtr w,
               FloatTexture2DPtr un) {

    unsigned int uw = u0->GetWidth();
    unsigned int uh = u0->GetHeight();
    FloatTexture2DPtr utp = FloatTexture2DPtr(new Texture2D<REAL>(uw,uh,1));

    REAL dx = 1, dy = 1;
    REAL beta = 0.00001;
    REAL my = 50000; //0.026691;
    
    REAL errorsum = 0;
    for (unsigned int i=1; i<uw-1; i++) {
        for (unsigned int k=1; k<uh-1; k++) {
            REAL gx = ((*un)(i+1,k) - (*un)(i-1,k)) / (2*dx);
            REAL gy = ((*un)(i,k+1) - (*un)(i,k-1)) / (2*dy);
            REAL gxx = ((*un)(i+1,k) - 2* (*un)(i,k) + (*un)(i-1,k))/(dx*dx);
            REAL gyy = ((*un)(i,k+1) - 2* (*un)(i,k) + (*un)(i,k-1))/(dx*dx);

            REAL gxy = ( (*un)(i+1,k+1) - (*un)(i-1,k+1) -
                         (*un)(i+1,k-1) + (*un)(i+1,k+1) )/(2*dx*dy);

            REAL ugx = 0;
            if ( (gx*( (*w)(i,k) - (*v0)(i,k))) > 0) {
                ugx = ((*un)(i,k) - (*un)(i-1,k) )/dx;
            } else if ((gx*( (*w)(i,k) - (*v0)(i,k))) < 0) {
                ugx = ((*un)(i+1,k) - (*un)(i,k) )/dx;
            }/* else {
                std::string error = "ugx is zero on index: " + 
                    Utils::Convert::ToString(i);
                error += ", " + Utils::Convert::ToString(k);
                throw Core::Exception(error);
                }*/
            
            REAL ugy = 0;
            if ( (gy*( (*w)(i,k) - (*v0)(i,k))) > 0) {
                ugy = ((*un)(i,k) - (*un)(i,k-1) )/dy;
            } else if ((gy*( (*w)(i,k) - (*v0)(i,k))) < 0) {
                ugy = ((*un)(i,k+1) - (*un)(i,k) )/dy;
            } /*else {
                std::string error = "ugy is zero on index: " + 
                    Utils::Convert::ToString(i);
                error += ", " + Utils::Convert::ToString(k);
                throw Core::Exception(error);
                }*/

            REAL sn = 0;
            if ( gx*gx + gy*gy > beta ) {
                sn = (gxx*gy*gy - 2*gxy*gx*gy + gyy*gx*gx)/(gx*gx+gy*gy);
            }

            REAL ugLength = 0;
            if ( ugx != 0 || ugy != 0)
                ugLength = -sqrt(ugx*ugx+ugy*ugy);

            REAL error = (-ugLength * ( (*w)(i,k)-(*v0)(i,k) ) +sn);
            errorsum += error;
            if (isnan(errorsum)) {
                logger.info <<"nan -exit at(" << i << "," << k 
                            << ")" << logger.end;
                exit(-1);
            }
            (*utp)(i,k) = error *my*dt + (*un)(i,k);
        }
    }

    logger.info << "errorsum: " << errorsum << logger.end;

    // solve boundary conditions
    for (unsigned int y=1; y<uh-1; y++) {
        (*utp)(0,y) = (*utp)(1,y);
        (*utp)(uw-1,y) = (*utp)(uw-2,y);
    }
    for (unsigned int x=0; x<uw; x++) {
        (*utp)(x,0) = (*utp)(x,1);
        (*utp)(x,uh-1) = (*utp)(x,uh-2);
    }
    
    return utp;
}

FloatTexture2DPtr ROFEquation(REAL dt,
                              FloatTexture2DPtr u0,
                              FloatTexture2DPtr v0) {
    
    FloatTexture2DPtr ut(u0->Clone());
    for (unsigned int i=0; i<40/*180*/; i++) {
        FloatTexture2DPtr Jt = J(dt, ut);
        FloatTexture2DPtr w = J(dt, FloatTexture2DPtr(Jt->Clone()));

        FloatTexture2DPtr utp = ROFIteration(dt,u0,v0,w,ut);

        logger.info << "ROF Equation (iteration " << i << "), " << logger.end;
    
        // swap
        FloatTexture2DPtr tmp = ut;
        ut = utp;
        utp = tmp;
    }
    return ut;
}

int main(int argc, char** argv) {
    // timer to mesure execution time
    Utils::Timer timer;
    timer.Start();

    // create a logger to std out
    StreamLogger* stdlog = new StreamLogger(&std::cout);
    Logger::AddLogger(stdlog);

    // prepare for resource loading
    DirectoryManager::AppendPath("./projects/DeConv/data/");
    ResourceManager<ITextureResource>::AddPlugin(new FreeImagePlugin());

    // load the RGBA gray scale input file
    string filename = "lena.png";
    ITexture2DPtr tex =
        ResourceManager<ITexture2D>::Create(filename);
    tex->Load();    

    // peel out one channel for gray scale from the RGBA texture.
    EmptyTextureResourcePtr ut_tex =
        EmptyTextureResource::CloneChannel(tex,0);
    //@todo: RGBAToGrayScale()

    // convert to floating point arrays
    FloatTexture2DPtr u0 = Utils::TexUtils::ToFloatTexture(ut_tex);

    // create a output canvas
    UCharTexture2DPtr output = Utils::TexUtils::ToUCharTexture(u0);
    TextureTool<unsigned char>::DumpTexture(Utils::TexUtils::ToRGBAfromLuminance(output), "u0.png");

    // run initiali heat equation
    FloatTexture2DPtr v0 = J(0.24, u0);
    output = Utils::TexUtils::ToUCharTexture(v0);
    TextureTool<unsigned char>::DumpTexture(Utils::TexUtils::ToRGBAfromLuminance(output), "v0.png");

    FloatTexture2DPtr result = ROFEquation(0.000001, u0, v0);
    string outputFile = "output.png";
    output = Utils::TexUtils::ToUCharTexture(result);
    TextureTool<unsigned char>::DumpTexture(Utils::TexUtils::ToRGBAfromLuminance(output), outputFile);

    logger.info << "execution time: " << timer.GetElapsedTime() << logger.end;
    return EXIT_SUCCESS;
}
