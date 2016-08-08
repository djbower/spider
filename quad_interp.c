/* adapted from:
       https://anisrahman.wordpress.com/2010/02/02/bilinear-interpolation/

original version is below this function */

int BilinearInterpolation_quad(float *in, float *vx, float *vy, int w, int h, float *out)
{       
    float fraction_x, fraction_y, one_minus_x, one_minus_y;
    int ceil_x, ceil_y, floor_x, floor_y;
 
    float pix[4];
 
    for (int y = 0; y < h ; ++y)
        for (int x = 0; x < w ; ++x)
        {
            if (vx[y*w+x]<0.0f || vx[y*w+x]>w-1 
             || vy[y*w+x]<0.0f || vy[y*w+x]>h-1) 
             { out[y*w+x] = 0.0f; continue; }
 
            floor_x = (int)floor(vx[y*w+x]);            
            floor_y = (int)floor(vy[y*w+x]);
 
            if (floor_x < 0) floor_x = 0;
            if (floor_y < 0) floor_y = 0;
 
            ceil_x = floor_x + 1;
            if (ceil_x >= w) ceil_x = floor_x;
 
            ceil_y = floor_y + 1;
            if (ceil_y >= h) ceil_y = floor_y;
             
            fraction_x = vx[y*w+x] - (float)floor_x;
            fraction_y = vy[y*w+x] - (float)floor_y;
             
            one_minus_x = 1.0 - fraction_x;
            one_minus_y = 1.0 - fraction_y;
 
            pix[0] = in[floor_y*w + floor_x];
            pix[1] = in[floor_y*w + ceil_x];
            pix[2] = in[ceil_y*w + floor_x];
            pix[3] = in[ceil_y*w + ceil_x];
 
            out[y*w + x] = 
one_minus_y*(one_minus_x*pix[0]+fraction_x*pix[1]) +
 fraction_y*(one_minus_x*pix[2]+fraction_x*pix[3]);
        }
 
        return 0;
}

/* original from:
       https://anisrahman.wordpress.com/2010/02/02/bilinear-interpolation/ */

int BilinearInterpolation(float *in, float *vx, float *vy, int w, int h, float *out)
{       
    float fraction_x, fraction_y, one_minus_x, one_minus_y;
    int ceil_x, ceil_y, floor_x, floor_y;
 
    float pix[4];
 
    for (int y = 0; y < h ; ++y)
        for (int x = 0; x < w ; ++x)
        {
            if (vx[y*w+x]<0.0f || vx[y*w+x]>w-1 
             || vy[y*w+x]<0.0f || vy[y*w+x]>h-1) 
             { out[y*w+x] = 0.0f; continue; }
 
            floor_x = (int)floor(vx[y*w+x]);            
            floor_y = (int)floor(vy[y*w+x]);
 
            if (floor_x < 0) floor_x = 0;
            if (floor_y < 0) floor_y = 0;
 
            ceil_x = floor_x + 1;
            if (ceil_x >= w) ceil_x = floor_x;
 
            ceil_y = floor_y + 1;
            if (ceil_y >= h) ceil_y = floor_y;
             
            fraction_x = vx[y*w+x] - (float)floor_x;
            fraction_y = vy[y*w+x] - (float)floor_y;
             
            one_minus_x = 1.0 - fraction_x;
            one_minus_y = 1.0 - fraction_y;
 
            pix[0] = in[floor_y*w + floor_x];
            pix[1] = in[floor_y*w + ceil_x];
            pix[2] = in[ceil_y*w + floor_x];
            pix[3] = in[ceil_y*w + ceil_x];
 
            out[y*w + x] = 
one_minus_y*(one_minus_x*pix[0]+fraction_x*pix[1]) +
 fraction_y*(one_minus_x*pix[2]+fraction_x*pix[3]);
        }
 
        return 0;
}
