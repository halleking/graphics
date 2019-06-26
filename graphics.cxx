#include <iostream>
#include <cstring>
#include <utility>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <math.h>
#define NORMALS

using std::cerr;
using std::endl;

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Triangle
{
  public:
      double          X[3];
      double          Y[3];
      double 		  Z[3];
      double    shading[3];
      double  colors[3][3];
      double normals[3][3];

};

class Screen
{
  public:
      unsigned char   *buffer;
      double	  *zbuffer;
      int width, height, npixels;

      void rasterize(int index, unsigned char RGBval)
      {
          this->buffer[index] = RGBval;
      }
};

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix          ViewTransform(Camera c);
    Matrix          CameraTransform(Camera c);
    Matrix          DeviceTransform(Screen s, Camera c);
};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

double Norm(double vector[3])
{
    double norm = sqrt((vector[0] * vector[0]) + (vector[1] * vector[1]) + (vector[2] * vector[2]));
    return norm;
}
double Dot(double v1[3], double v2[3])
{
    return ((v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2]));
}

double *CrossProduct(double v1[3], double v2[3])
{
	double *cp = new double[3];
	cp[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
	cp[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
	cp[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
	return cp;
}

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

Matrix Camera::CameraTransform(Camera camera)
{
    //CALCULATE CAMERA FRAME
    double *o = new double[3];
    double *w = new double[3];
    double *u = new double[3];
    double *v = new double[3];
    double *t = new double[3];
    o = camera.position;
    w[0] = o[0] - camera.focus[0];
    w[1] = o[1] - camera.focus[1];
    w[2] = o[2] - camera.focus[2];
    u = CrossProduct(camera.up, w);
    v = CrossProduct(w, u);
    double nw = Norm(w);
    double nu = Norm(u);
    double nv = Norm(v);
    for (int i = 0; i < 3; i++)
    {
    	w[i] = w[i] / nw;
    	u[i] = u[i] / nu;
    	v[i] = v[i] / nv;
    }
   
    t[0] = -o[0];
    t[1] = -o[1];
    t[2] = -o[2];

    Matrix cf;
    cf.A[0][0] = u[0];
    cf.A[0][1] = v[0];
    cf.A[0][2] = w[0];
    cf.A[0][3] = 0;
    cf.A[1][0] = u[1];
    cf.A[1][1] = v[1];
    cf.A[1][2] = w[1];
    cf.A[1][3] = 0;
    cf.A[2][0] = u[2];
    cf.A[2][1] = v[2];
    cf.A[2][2] = w[2];
    cf.A[2][3] = 0;
    cf.A[3][0] = Dot(u,t);
    cf.A[3][1] = Dot(v,t);
    cf.A[3][2] = Dot(w,t);
    cf.A[3][3] = 1;
    return cf;
};

Matrix Camera::ViewTransform(Camera camera)
{
    double a = camera.angle;
    double n = camera.near;
   	double f = camera.far;

    Matrix vf;
    vf.A[0][0] = 1 / tan(a/2);
    vf.A[0][1] = 0;
    vf.A[0][2] = 0;
    vf.A[0][3] = 0;
    vf.A[1][0] = 0;
    vf.A[1][1] = 1 / tan(a/2);
    vf.A[1][2] = 0;
    vf.A[1][3] = 0;
    vf.A[2][0] = 0;
    vf.A[2][1] = 0;
    vf.A[2][2] = (f+n) / (f-n);
    vf.A[2][3] = -1;
    vf.A[3][0] = 0;
    vf.A[3][1] = 0;
    vf.A[3][2] = (2*f*n) / (f-n);
    vf.A[3][3] = 0;
    return vf;
}

Matrix Camera::DeviceTransform(Screen s, Camera c)
{

	double n = s.width;
	double m = s.height;
    Matrix dt;
    dt.A[0][0] = n/2;
    dt.A[0][1] = 0;
    dt.A[0][2] = 0;
    dt.A[0][3] = 0;
    dt.A[1][0] = 0;
    dt.A[1][1] = m/2;
    dt.A[1][2] = 0;
    dt.A[1][3] = 0;
    dt.A[2][0] = 0;
    dt.A[2][1] = 0;
    dt.A[2][2] = 1;
    dt.A[2][3] = 0;
    dt.A[3][0] = n/2;
    dt.A[3][1] = m/2;
    dt.A[3][2] = 0;
    dt.A[3][3] = 1;
    return dt;
}

void RasterizeBottom(Screen s, Triangle t)
{
    int TL, TR, B;
    //X values
    double x1 = t.X[0];
    double x2 = t.X[1];
    double x3 = t.X[2];
    //Y values
    double y1 = t.Y[0];
    double y2 = t.Y[1];
    double y3 = t.Y[2];

    double minY = std::min(std::min(y1,y2),y3);
    double maxY = std::max(std::max(y1,y2),y3);

    if (y1 == maxY)
    {
        if (y2 == maxY)
        {
            if (x1 < x2)
            {
                TL = 0; 
                TR = 1;
                B = 2;
            }
           else
           {
                TL = 1; 
                TR = 0;
                B = 2;
            }
        }
        else if (y3 == maxY)
        {
            if (x1 < x3){
                TL = 0;
                TR = 2;
                B = 1;
            }
            else
            {
                TL = 2;
                TR = 0;
                B = 1;
            }
        }
    }
    else if (y2 == maxY)
    {
        if (x2 < x3)
        {
            TL = 1; 
            TR = 2;
            B = 0;
        }
        else
        {
            TL = 2; 
            TR = 1;
            B = 0;
        }
    }
    //DETERMINE ROWS OF PIXELS TRIANGLES CAN INTERSECT
    double rowMin = ceil_441(minY);
    double rowMax = floor_441(maxY);
   
    double leftX = t.X[TL];
    double leftY = t.Y[TL];
    double rightX = t.X[TR];
    double rightY = t.Y[TR];
    double botX = t.X[B];
    double botY = t.Y[B];

    for (int r = rowMin; r <= rowMax; r++)
    {
        //CHECK BOUNDS OF ROWS
        if (r < 0 || r > s.height)
            continue;

        //FIND LEFT END
        double leftEnd;
        double leftSlope;
        double intercept;
        double ldx = (botX - leftX);
        double ldy = (botY - leftY);
        if (ldx == 0)
        {
            leftSlope = 0;
            leftEnd = leftX;
        }
        else
        {
            leftSlope = ldy / ldx;
            intercept = leftY - (leftSlope * leftX);
            leftEnd = (r - intercept) / leftSlope;
        }
        
        //FIND RIGHT END
        double rightEnd;
        double rightSlope;
        double intercept1;
        double rdx = (botX - rightX);
        double rdy = (botY - rightY);
        if (rdx == 0)
        {
            rightSlope = 0;
            rightEnd = rightX;
        }
        else
        {
            rightSlope = rdy / rdx;
            intercept1 = rightY - (rightSlope * rightX);
            rightEnd = (r - intercept1) / rightSlope;
        }

        //LERP TO FIND LEFT AND RIGHT ENDS
        double yt = (r - botY) / (leftY - botY);

        double LRed = t.colors[B][0] + yt * (t.colors[TL][0] - t.colors[B][0]);
        double LGreen = t.colors[B][1] + yt * (t.colors[TL][1] - t.colors[B][1]);
        double LBlue = t.colors[B][2] + yt * (t.colors[TL][2] - t.colors[B][2]);
        double LZ = t.Z[B] + yt * (t.Z[TL] - t.Z[B]);
        double LS = t.shading[B] + yt * (t.shading[TL] - t.shading[B]);

        double RRed = t.colors[B][0] + yt * (t.colors[TR][0] - t.colors[B][0]);
        double RGreen = t.colors[B][1] + yt * (t.colors[TR][1] - t.colors[B][1]);
        double RBlue = t.colors[B][2] + yt * (t.colors[TR][2] - t.colors[B][2]);
        double RZ = t.Z[B] + yt * (t.Z[TR] - t.Z[B]);
        double RS = t.shading[B] + yt * (t.shading[TR] - t.shading[B]);

        //SCANLINE COLUMNS OF PIXELS THAT CAN INTERSECT 
        double ceilLeft = ceil_441(leftEnd);
        double floorRight = floor_441(rightEnd);     

        for (int c = ceilLeft; c <= floorRight; c++)
        {
            //CHECK BOUNDS OF COLUMNS 
            if (c < 0 || c >= s.width) 
              continue;

          	//GET INDEX OF PIXEL
            int index = 3 * (r*s.width+c);

            //CHECK THAT ALL PIXELS ARE ON SCREEN
            if (index < 0 || index+2 > 3*s.npixels)
                continue;

            //LERP BETWEEN LEFT AND RIGHT ENDS
            double xt = (c - leftEnd) / (rightEnd - leftEnd);
            double Red = LRed + xt * (RRed - LRed);
            double Green = LGreen + xt * (RGreen - LGreen);
            double Blue = LBlue + xt * (RBlue - LBlue);
            double Z = LZ + xt * (RZ - LZ);
            double S = LS + xt * (RS - LS);

            if (Z > s.zbuffer[r*s.width+c])
            {
            	s.zbuffer[r*s.width+c] = Z;
            	s.rasterize(index, ceil_441(255*fmin(1,Red * S)));
            	s.rasterize(index+1, ceil_441(255*fmin(1,Green * S)));
            	s.rasterize(index+2, ceil_441(255*fmin(1,Blue * S)));
            }
        }
    }
}

void RasterizeTop(Screen s, Triangle t)
{
    int BL, BR, T;
    double ceilLeft, floorRight;

    //X values
    double x1 = t.X[0];
    double x2 = t.X[1];
    double x3 = t.X[2];

    //Y values
    double y1 = t.Y[0];
    double y2 = t.Y[1];
    double y3 = t.Y[2];

    double minY = std::min(std::min(y1,y2),y3);
    double maxY = std::max(std::max(y1,y2),y3);

    if (y1 == minY)
    {
        if (y2 == minY)
        {
            if (x1 < x2)
            {
                BL = 0; 
                BR = 1;
                T = 2;
            }
           else
           {
                BL = 1; 
                BR = 0;
                T = 2;
            }
        }
        else if (y3 == minY)
        {
            if (x1 < x3){
                BL = 0;
                BR = 2;
                T = 1;
            }
            else
            {
                BL = 2;
                BR = 0;
                T = 1;
            }
        }
    }
    else if (y2 == minY)
    {
        if (x2 < x3)
        {
            BL = 1; 
            BR = 2;
            T = 0;
        }
        else
        {
            BL = 2; 
            BR = 1;
            T = 0;
        }
    }
    //DETERMINE ROWS OF PIXELS TRIANGLES CAN INTERSECT
    double rowMin = ceil_441(minY);
    double rowMax = floor_441(maxY);

    double leftX = t.X[BL];
    double leftY = t.Y[BL];
    double rightX = t.X[BR];
    double rightY = t.Y[BR];
    double topX = t.X[T];
    double topY = t.Y[T];

    for (int r = rowMin; r <= rowMax; r++)
    {
        //CHECK BOUNDS OF ROWS
        if (r < 0 || r > s.height)
            continue;

        //FIND LEFT END
        double leftEnd;
        double leftSlope;
        double intercept;
        double ldx = (topX - leftX);
        double ldy = (topY - leftY);
        if (ldx == 0)
        {
            leftSlope = 0;
            leftEnd = leftX;
        }
        else
        {
            leftSlope = ldy / ldx;
            intercept = leftY - (leftSlope * leftX);
            leftEnd = (r - intercept) / leftSlope;
        }
  
        //FIND RIGHT END
        double rightEnd;
        double rightSlope;
        double intercept1;
        double rdx = (topX - rightX);
        double rdy = (topY - rightY);
        if (rdx == 0)
        {
            rightSlope = 0;
            rightEnd = rightX;
        }

        else
        {
            rightSlope = rdy / rdx;
            intercept1 = rightY - (rightSlope * rightX);
            rightEnd = (r - intercept1) / rightSlope;
        }
        //LERP TO FIND LEFT AND RIGHT ENDS
        double yt = (r - topY) / (leftY - topY);

        double LRed = t.colors[T][0] + yt * (t.colors[BL][0] - t.colors[T][0]);
        double LGreen = t.colors[T][1] + yt * (t.colors[BL][1] - t.colors[T][1]);
        double LBlue = t.colors[T][2] + yt * (t.colors[BL][2] - t.colors[T][2]);
        double LZ = t.Z[T] + yt * (t.Z[BL] - t.Z[T]);
        double LS = t.shading[T] + yt * (t.shading[BL] - t.shading[T]);

        double RRed = t.colors[T][0] + yt * (t.colors[BR][0] - t.colors[T][0]);
        double RGreen = t.colors[T][1] + yt * (t.colors[BR][1] - t.colors[T][1]);
        double RBlue = t.colors[T][2] + yt * (t.colors[BR][2] - t.colors[T][2]);
        double RZ = t.Z[T] + yt * (t.Z[BR] - t.Z[T]);
        double RS = t.shading[T] + yt * (t.shading[BR] - t.shading[T]);

        //SCANLINE COLUMNS OF PIXELS THAT CAN INTERSECT 
        ceilLeft = ceil_441(leftEnd);
        floorRight = floor_441(rightEnd); 
        
        for (int c = ceilLeft; c <= floorRight; c++)
        {
            //CHECK BOUNDS OF COLUMNS  
            if (c < 0 || c >= s.width) 
                continue;

            //FIND INDEX OF PIXEL
            int index = 3 * (r*s.width+c);

            //CHECK THAT ALL PIXELS ARE ON SCREEN
            if (index < 0 || index+2 > 3*s.npixels)           
                continue;

            //LERP BETWEEN LEFT AND RIGHT ENDS
            double xt = (c - leftEnd) / (rightEnd - leftEnd);
            double Red = LRed + xt * (RRed - LRed);
            double Green = LGreen + xt * (RGreen - LGreen);
            double Blue = LBlue + xt * (RBlue - LBlue);
            double Z = LZ + xt * (RZ - LZ);
            double S = LS + xt * (RS - LS);

            if (Z > s.zbuffer[r*s.width+c])
            {
            	s.zbuffer[r*s.width+c] = Z;
            	s.rasterize(index, ceil_441(255*fmin(1,Red * S)));
            	s.rasterize(index+1, ceil_441(255*fmin(1,Green * S)));
            	s.rasterize(index+2, ceil_441(255*fmin(1,Blue * S)));
            }
        }
    }
}

void Render(std::vector<Triangle> triangles, Screen screen)
{
	//LOOP THROUGH ALL TRIANGLES
   for (int i = 0; i < triangles.size(); i++)
   {
        //X values
        double x1 = triangles[i].X[0];
        double x2 = triangles[i].X[1];
        double x3 = triangles[i].X[2];

        //Y values
        double y1 = triangles[i].Y[0];
        double y2 = triangles[i].Y[1];
        double y3 = triangles[i].Y[2];

        int mid, top, bot;

        //ARBITRARY TRIANGLE
        if (y1 != y2 && y1 != y3 && y2 != y3)
        { 
            if (y1 > y2 && y1 > y3)
            {
                top = 0;
                if (y2 > y3)
                {
                    mid = 1;
                    bot = 2;
                }
                else
                {
                    mid = 2;
                    bot = 1;
                }
            }
            else if (y2 > y1 && y2 > y3)
            {
                top = 1;
                if (y1 > y3)
                {
                    mid = 0; 
                    bot = 2;
                }
                else
                {
                    mid = 2; 
                    bot = 0;
                }
            }
            else if (y3 > y1 && y3 > y2)
            {
                top = 2; 
                if (y1 > y2)
                {
                    mid = 0;
                    bot = 1;
                }
                else
                {
                    mid = 1;
                    bot = 0;
                }
            }
            
            // SLOPE BETWEEN TOP AND BOTTOM POINTS
            double ldx = (triangles[i].X[bot] - triangles[i].X[top]);
            double ldy = (triangles[i].Y[bot] - triangles[i].Y[top]);
            double slope;
            double intercept;
            double x;

            if (ldx == 0)
            {
              slope = 0;
              x = triangles[i].X[top];
            }
            else
            {
              slope = ldy / ldx;
              intercept = triangles[i].Y[top] - (slope * triangles[i].X[top]);
              x = (triangles[i].Y[mid] - intercept) / slope;
            }

            //LERP for mid Z 
            double t = (triangles[i].Y[mid] - triangles[i].Y[bot]) / (triangles[i].Y[top] - triangles[i].Y[bot]);
            double newZ = triangles[i].Z[bot] + t * (triangles[i].Z[top] - triangles[i].Z[bot]);

            //LERP for mid color
            double newR = triangles[i].colors[bot][0] + t * (triangles[i].colors[top][0] - triangles[i].colors[bot][0]);
            double newG = triangles[i].colors[bot][1] + t * (triangles[i].colors[top][1] - triangles[i].colors[bot][1]);
            double newB = triangles[i].colors[bot][2] + t * (triangles[i].colors[top][2] - triangles[i].colors[bot][2]);

            //LERP for mid shading
            double newS = triangles[i].shading[bot] + t * (triangles[i].shading[top] - triangles[i].shading[bot]);

            //MAKE TOP TRIANGLE
            Triangle TOP; 
            double topX[3] = {x, triangles[i].X[top], triangles[i].X[mid]};
            double topY[3] = {triangles[i].Y[mid], triangles[i].Y[top], triangles[i].Y[mid]};
            double topZ[3] = {newZ, triangles[i].Z[top], triangles[i].Z[mid]};
            double topS[3] = {newS, triangles[i].shading[top], triangles[i].shading[mid]};
            double topColor[3][3] = {{newR, newG, newB}, {triangles[i].colors[top][0], triangles[i].colors[top][1], triangles[i].colors[top][2]}, {triangles[i].colors[mid][0], triangles[i].colors[mid][1], triangles[i].colors[mid][2]}};
            memcpy(TOP.X, topX, sizeof(double)*3);
            memcpy(TOP.Y, topY, sizeof(double)*3);
            memcpy(TOP.Z, topZ, sizeof(double)*3);
            memcpy(TOP.shading,topS, sizeof(double)*3);
            memcpy(TOP.colors, topColor, sizeof(double[3])*3); 

            //MAKE BOTTOM TRIANGLE
            Triangle BOTTOM; 
            double botX[3] = {x, triangles[i].X[bot], triangles[i].X[mid]};
            double botY[3] = {triangles[i].Y[mid], triangles[i].Y[bot], triangles[i].Y[mid]};
            double botZ[3] = {newZ, triangles[i].Z[bot], triangles[i].Z[mid]};
            double botS[3] = {newS, triangles[i].shading[bot], triangles[i].shading[mid]};
            double botColor[3][3] = {{newR, newG, newB}, {triangles[i].colors[bot][0], triangles[i].colors[bot][1], triangles[i].colors[bot][2]}, {triangles[i].colors[mid][0], triangles[i].colors[mid][1], triangles[i].colors[mid][2]}};
            memcpy(BOTTOM.X, botX, sizeof(double)*3);
            memcpy(BOTTOM.Y, botY, sizeof(double)*3);
            memcpy(BOTTOM.Z, botZ, sizeof(double)*3);
            memcpy(BOTTOM.shading,botS, sizeof(double)*3);
            memcpy(BOTTOM.colors, botColor, sizeof(double[3])*3);

            //RASTERIZE SPLIT TRIANGLES
            RasterizeBottom(screen, BOTTOM);
           	RasterizeTop(screen, TOP);

      }
      //UP OR DOWN TRIANGLES
      else if (y1 == y2 || y1 == y3 || y2 == y3)
      {
      		if (y1 == y2)
      		{
      			if (y3 < y1)
      				RasterizeBottom(screen, triangles[i]);
      			else
      				RasterizeTop(screen, triangles[i]);
      		}
      		else if (y1 == y3)
      		{
      			if (y2 < y1)
      				RasterizeBottom(screen, triangles[i]);
      			else
      				RasterizeTop(screen, triangles[i]);
      		}
      		else if (y2 == y3)
      		{
      			if (y1 < y2)
      				RasterizeBottom(screen, triangles[i]);
      			else
      				RasterizeTop(screen, triangles[i]);
      		}
      }
   }
}


double *CalculateShading(Triangle t, Camera c)
{
    double *shading = new double[3];
    for (int i = 0; i < 3; i++)
    {
        double diffuse, specular;
        double LdotN = Dot(lp.lightDir, t.normals[i]);
        double *r = new double[3];
        for (int j = 0; j < 3; j++)
        {
             r[j] = (2*LdotN * t.normals[i][j]) - lp.lightDir[j];
        }
        //CALCULATE VIEWDIR
        //ViewingDirection = CameraPosition - VertexPosition
        double *viewDir = new double[3]; 
        viewDir[0] = c.position[0] - t.X[i];
        viewDir[1] = c.position[1] - t.Y[i];
        viewDir[2] = c.position[2] - t.Z[i];

        //NORMALIZE VECTORS
        double lpn = Norm(lp.lightDir);
        double nn = Norm(t.normals[i]);
        double rn = Norm(r);
        double vn = Norm(viewDir);
        for (int k = 0; k < 3; k++)
        {
            lp.lightDir[k] = lp.lightDir[k] / lpn;
            t.normals[i][k] = t.normals[i][k] / nn;
            r[k] = r[k] / rn;
            viewDir[k] = viewDir[k] / vn;
        }

        //CALCULATE SPECULAR SHADING
        double RdotV = Dot(r, viewDir);
        if (RdotV < 0){
            specular = 0;
        }
        else{
            double powg = pow(RdotV, lp.alpha);
            specular = std::max(0.0, lp.Ks * powg);
        }

        //CALCULATE DIFFUSE SHADING (2-sided)
        diffuse = lp.Kd * abs(LdotN);

        //SET SHADING AMOUNT FOR VERTEX
        double phong = lp.Ka + diffuse + specular;
        shading[i] = phong;
    } 

    return shading;
}


void Transform(Screen screen, Camera camera, std::vector<Triangle> triangles)
{
    Matrix ct = camera.CameraTransform(camera);
    Matrix vt = camera.ViewTransform(camera);
    Matrix dt = camera.DeviceTransform(screen, camera);
    Matrix M = Matrix::ComposeMatrices(Matrix::ComposeMatrices(ct, vt), dt);

    //MULTIPLY VERTICES BY M TO GET TRANSFORM
    std::vector<Triangle> transformed;
    transformed = triangles;
    for (int i = 0; i < transformed.size(); i++)
    {
        double *shading = CalculateShading(transformed[i], camera);

        for (int j = 0; j < 3; j++)
        {
            transformed[i].shading[j] = shading[j];
        }

    	double *v0 = new double[4];
    	v0[0] = transformed[i].X[0];
    	v0[1] = transformed[i].Y[0];
    	v0[2] = transformed[i].Z[0];
    	v0[3] = 1;
    	double *v1 = new double[4];
    	v1[0] = transformed[i].X[1];
    	v1[1] = transformed[i].Y[1];
    	v1[2] = transformed[i].Z[1];
    	v1[3] = 1;
    	double *v2 = new double[4];
    	v2[0] = transformed[i].X[2];
    	v2[1] = transformed[i].Y[2];
    	v2[2] = transformed[i].Z[2];
    	v2[3] = 1;
    	
    	double *newv = new double[4];

    	M.TransformPoint(v0, newv);
    	transformed[i].X[0] = newv[0]/ newv[3];
    	transformed[i].Y[0] = newv[1]/ newv[3];
    	transformed[i].Z[0] = newv[2]/ newv[3];
    	
    	M.TransformPoint(v1, newv);
    	transformed[i].X[1] = newv[0]/ newv[3];
    	transformed[i].Y[1] = newv[1]/ newv[3];
    	transformed[i].Z[1] = newv[2]/ newv[3];

    	M.TransformPoint(v2, newv);
    	transformed[i].X[2] = newv[0]/ newv[3];
    	transformed[i].Y[2] = newv[1]/ newv[3];
    	transformed[i].Z[2] = newv[2]/ newv[3];
    }
    //APPLY RASTERIZATION/ZBUFFER
    Render(transformed, screen);
}

int main()
{
   	//GET TRIANGLES FROM FILE
    std::vector<Triangle> TRIANGLES = GetTriangles();

   //LOOP THROUGH DIFFERENT CAMERA POSITIONS
   for (int i = 0; i < 1000; i++)
   {
   		//INITIALIZE SCREEN
   		vtkImageData *image = NewImage(1000, 1000);
   		int npixels = 1000*1000;
   		unsigned char *buffer = 
     		(unsigned char *) image->GetScalarPointer(0,0,0);
   		double *zbuffer = new double[npixels];
   		//INITIALIZE BUFFERS
   		for (int i = 0 ; i < npixels*3 ; i++)
       		buffer[i] = 0;
   		for (int j = 0; j < npixels; j++)
   				zbuffer[j] = -1;
   		Screen screen;
   		screen.buffer = buffer;
   		screen.zbuffer = zbuffer;
   		screen.npixels = npixels;
   		screen.width = 1000;
   		screen.height = 1000;

   		//GET CAMERA POSITION
   		Camera c = GetCamera(i, 1000);

   		//TRANSFORM TO DEVICE SPACE & RENDER
   		Transform(screen, c, TRIANGLES);

   		//SAVE IMAGE
        char filename[32];
        sprintf(filename, "frame%03d", i);
        WriteImage(image, filename);
   }
}
