#ifndef GRADIENT_H
#define GRADIENT_H

/*----------------------------------------------------------------------------*/
/*--------------------------------- Gradient ---------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Computes the direction of the level line of 'in' at each point2i.

    The result is:
    - an image_double with the angle at each pixel, or NOTDEF if not defined.
    - the image_double 'modgrad' (a point2ier is passed as argument)
      with the gradient magnitude at each point2i.
    - a list of pixels 'list_p' roughly ordered by decreasing
      gradient magnitude. (The order is made by classifying point2is
      into bins by gradient magnitude. The parameters 'n_bins' and
      'max_grad' specify the number of bins and the gradient modulus
      at the highest bin. The pixels in the list would be in
      decreasing gradient magnitude, up to a precision of the size of
      the bins.)
    - a point2ier 'mem_p' to the memory used by 'list_p' to be able to
      free the memory when it is not used anymore.
 */
//返回一张梯度角度顺时针旋转90°后的align角度图angles，如果梯度角度是(gx,gy)->(-gy,gx)，
//和梯度的模的图modgrad,然后按照n_bins进行伪排序返回链表的头指针list_p,里面存的是坐标
//(Returns an align angle map angles after the gradient angle is rotated 90° clockwise, 
//if the gradient angle is (gx,gy)->(-gy,gx), and the map modgrad of the gradient's modulus, 
//and then returns according to the pseudo-sort of n_bins The head pointer list_p of the linked list, 
//which stores the coordinates)
static image_double ll_angle( image_double in, double threshold,
                              struct coorlist ** list_p,
                              image_double * modgrad, unsigned int n_bins )
{
  image_double g;
  unsigned int n,p,x,y,adr,i;
  double com1,com2,gx,gy,norm,norm2;
  /* the rest of the variables are used for pseudo-ordering
     the gradient magnitude values */
  int list_count = 0;
  //struct coorlist * list;
  struct coorlist *temp;
  struct coorlist ** range_l_s; /* array of point2iers to start of bin list,表示1024个bin的头指针的指针数组 */
  struct coorlist ** range_l_e; /* array of point2iers to end of bin list，表示1024个bin的尾指针的指针数组*/
  struct coorlist * start;
  struct coorlist * end;
  double max_grad = 0.0;

  /* check parameters */
  if( in == NULL || in->data == NULL || in->xsize == 0 || in->ysize == 0 )
    error("ll_angle: invalid image.");
  if( threshold < 0.0 ) error("ll_angle: 'threshold' must be positive.");
  if( list_p == NULL ) error("ll_angle: NULL point2ier 'list_p'.");
 // if( mem_p == NULL ) error("ll_angle: NULL point2ier 'mem_p'.");
  if( modgrad == NULL ) error("ll_angle: NULL point2ier 'modgrad'.");
  if( n_bins == 0 ) error("ll_angle: 'n_bins' must be positive.");

  /* image size shortcuts */
  n = in->ysize;
  p = in->xsize;

  /* allocate output image */
  g = new_image_double(in->xsize,in->ysize);

  /* get memory for the image of gradient modulus */
  *modgrad = new_image_double(in->xsize,in->ysize);

  /* get memory for "ordered" list of pixels */
  //list = (struct coorlist *) calloc( (size_t) (n*p), sizeof(struct coorlist) );
  //*mem_p = (void *) list;
  range_l_s = (struct coorlist **) calloc( (size_t) n_bins,
                                           sizeof(struct coorlist *) );
  range_l_e = (struct coorlist **) calloc( (size_t) n_bins,
                                           sizeof(struct coorlist *) );
 // if( list == NULL || range_l_s == NULL || range_l_e == NULL )
  if( range_l_s == NULL || range_l_e == NULL )
    error("not enough memory.");
  for(i=0;i<n_bins;i++) range_l_s[i] = range_l_e[i] = NULL;

  /* 'undefined' on the down and right boundaries */
  for(x=0;x<p;x++) g->data[(n-1)*p+x] = NOTDEF;// p = in->xsize
  for(y=0;y<n;y++) g->data[p*y+p-1]   = NOTDEF;// n = in->ysize;

  /* compute gradient on the remaining pixels */
  for(x=0;x<p-1;x++)
    for(y=0;y<n-1;y++)
      {
        adr = y*p+x;

        /*
           Norm 2 computation using 2x2 pixel window:
             A B
             C D
           and
             com1 = D-A,  com2 = B-C.
           Then
             gx = B+D - (A+C)   horizontal difference
             gy = C+D - (A+B)   vertical difference
           com1 and com2 are just to avoid 2 additions.
         */
        com1 = in->data[adr+p+1] - in->data[adr];
        com2 = in->data[adr+1]   - in->data[adr+p];

        gx = com1+com2; /* gradient x component */
        gy = com1-com2; /* gradient y component */
        norm2 = gx*gx+gy*gy;
        norm = sqrt( norm2 / 4.0 ); /* gradient norm */

        (*modgrad)->data[adr] = norm; /* store gradient norm */

        if( norm <= threshold ) /* norm too small, gradient no defined */
          g->data[adr] = NOTDEF; /* gradient angle not defined */
        else
          {
            /* gradient angle computation */
            g->data[adr] = atan2(gx,-gy);

            /* look for the maximum of the gradient */
            if( norm > max_grad ) max_grad = norm;
          }
      }

  /* compute histogram of gradient values */
  for(x=0;x<p-1;x++)
    for(y=0;y<n-1;y++)
      {
		temp = new coorlist();
		if(temp == NULL)
		{
			printf("not enough memory");
			system("pause");
		}
        norm = (*modgrad)->data[y*p+x];
        /* store the point2i in the right bin according to its norm */
        i = (unsigned int) (norm * (double) n_bins / max_grad);
        if( i >= n_bins ) i = n_bins-1;
        if( range_l_e[i] == NULL )
          range_l_s[i] = range_l_e[i] = temp;//记录第i个区域的头指针到range_l_s[i] (Record the head pointer of the i-th area to range_l_s[i])
        else
          {
            range_l_e[i]->next = temp;//第i个区域由尾指针range_l_e[i]完成勾链 (The i-th area is linked by the tail pointer range_l_e[i])
            range_l_e[i] = temp;
          }
        range_l_e[i]->x = (int) x;//将坐标(x,y)记录到第i个分区 (Record the coordinates (x,y) to the ith partition)
        range_l_e[i]->y = (int) y;
        range_l_e[i]->next = NULL;
      }

  /* Make the list of pixels (almost) ordered by norm value.
     It starts by the larger bin, so the list starts by the
     pixels with the highest gradient value. Pixels would be ordered
     by norm value, up to a precision given by max_grad/n_bins.
   */
  for(i=n_bins-1; i>0 && range_l_s[i]==NULL; i--);//找到第一个不为空的分区bin (Find the first partition bin that is not empty)
  start = range_l_s[i];
  end = range_l_e[i];
  if( start != NULL )
    while(i>0)
      {
        --i;
        if( range_l_s[i] != NULL )
          {
            end->next = range_l_s[i];
            end = range_l_e[i];
          }
      }
  *list_p = start;
 // *mem_p  = start;
  /* free memory */
  free( (void *) range_l_s );
  free( (void *) range_l_e );

  return g;
}
/*----------------------------------------------------------------------------*/
/** Is point2i (x,y) aligned to angle theta, up to precision 'prec'?
 */
static int isaligned( int x, int y, image_double angles, double theta,
                      double prec )
{
  double a;

  /* check parameters */
  if( angles == NULL || angles->data == NULL )
    error("isaligned: invalid image 'angles'.");
  if( x < 0 || y < 0 || x >= (int) angles->xsize || y >= (int) angles->ysize )
    error("isaligned: (x,y) out of the image.");
  if( prec < 0.0 ) error("isaligned: 'prec' must be positive.");

  /* angle at pixel (x,y) */
  a = angles->data[ x + y * angles->xsize ];

  /* pixels whose level-line angle is not defined
     are considered as NON-aligned */
  if( a == NOTDEF ) return FALSE;  /* there is no need to call the function
                                      'double_equal' here because there is
                                      no risk of problems related to the
                                      comparison doubles, we are only
                                      interested in the exact NOTDEF value */

  /* it is assumed that 'theta' and 'a' are in the range [-pi,pi] */
  theta -= a;
  if( theta < 0.0 ) theta = -theta;
  if( theta > M_3_2_PI )
    {
	//--------------------------------------
	//origin code
     /* theta -= M_2__PI;
      if( theta < 0.0 ) theta = -theta;*/
	//--------------------------------------
	  //-------------------------------------
	  //mycode
	  theta = M_2__PI-theta;
	  if(theta < 0.0) 
		 theta = -theta; 
	  //--------------------------------------
    }

  return theta <= prec;
}

//计算梯度，返回模和角度，同时模值太小的像素点直接抑制掉，赋值为NOTDEF
// (Calculate the gradient, return the modulus and angle, and directly suppress the pixels whose modulus value is too small, and assign the value to NOTDEF)
//mod、angles为了传值，是二级指针
// (mod and angles are secondary pointers in order to pass values)
void calculateGradient( const double * img_in, unsigned int imgx, unsigned int imgy,image_double * mod, image_double * angles)
{
	if(img_in == NULL || imgx == 0 || imgy == 0)
		error("calculateGradient error!");
	(*mod) = new_image_double(imgx,imgy);
	(*angles) = new_image_double(imgx,imgy);
	double threshold = 2/sin(22.5/180*M_PI);
	unsigned int x,y,adr;
	double com1,com2;
	double gx,gy;
	double norm,norm_square;
	double sum = 0;

	//double max_grad = 0.0;
	//边界初始为NOTDEF (Boundary is initially NOTDEF)
	for ( x = 0; x<imgx; x++) 
	{
		//(*angles)->data[x]=NOTDEF;
		(*angles)->data[(imgy-1)*imgx+x]=NOTDEF;
		//(*mod)->data[x]=NOTDEF;
		(*mod)->data[(imgy-1)*imgx+x]=NOTDEF;
	}
	for ( y = 0; y<imgy; y++) 
	{
		//(*angles)->data[y*imgx] = NOTDEF;
		(*angles)->data[y*imgx+imgx-1] = NOTDEF;
		//(*mod)->data[y*imgx] = NOTDEF;
		(*mod)->data[y*imgx+imgx-1] = NOTDEF;
	}
	 /* compute gradient on the remaining pixels */
	for(x=0;x<imgx-1;x++)
		for(y=0;y<imgy-1;y++)
		{
			adr = y*imgx+x;
		  /*
		     Norm 2 computation using 2x2 pixel window:
		       A B
		       C D
		     and
		       com1 = D-A,  com2 = B-C.
		     Then
		       gx = B+D - (A+C)   horizontal difference
		       gy = C+D - (A+B)   vertical difference
		     com1 and com2 are just to avoid 2 additions.
		   */
		  com1 = img_in[adr+imgx+1] - img_in[adr];
		  com2 = img_in[adr+1]   - img_in[adr+imgx];

		  gx = com1+com2; /* gradient x component */
		  gy = com1-com2; /* gradient y component */
		  norm_square = gx*gx+gy*gy;

		  norm = sqrt( norm_square / 4.0 ); /* gradient norm */

		  (*mod)->data[adr] = norm; /* store gradient norm */

		  if( norm <= threshold ) /* norm too small, gradient no defined */
		  {
		    (*angles)->data[adr] = NOTDEF; /* gradient angle not defined */
			(*mod)->data[adr] = NOTDEF;
		  }
		  else
		    {
		      /* gradient angle computation */
		      (*angles)->data[adr] = atan2(gx,-gy);
		    }
		}
}

void calculateGradient2( const double * img_in, unsigned int imgx, unsigned int imgy, image_double * angles)
{
	if(img_in == NULL || imgx == 0 || imgy == 0)
		error("calculateGradient error!");
	image_double mod = new_image_double(imgx,imgy);
	(*angles) = new_image_double(imgx,imgy);
	unsigned int x,y,adr;
	double com1,com2;
	double gx,gy;
	double norm,norm_square;
	double threshold;
	double sum = 0;
	double value;  
	//double max_grad = 0.0;
	//边界初始为NOTDEF (Boundary is initially NOTDEF)
	for ( x = 0; x<imgx; x++) 
	{
		(*angles)->data[x]=NOTDEF;
		(*angles)->data[(imgy-1)*imgx+x]=NOTDEF;
		(mod)->data[x]=NOTDEF;
		(mod)->data[(imgy-1)*imgx+x]=NOTDEF;
	}
	for ( y = 0; y<imgy; y++) 
	{
		(*angles)->data[y*imgx] = NOTDEF;
		(*angles)->data[y*imgx+imgx-1] = NOTDEF;
		(mod)->data[y*imgx] = NOTDEF;
		(mod)->data[y*imgx+imgx-1] = NOTDEF;
	}
	 /* compute gradient on the remaining pixels */
	for(x=1;x<imgx-1;x++)
		for(y=1;y<imgy-1;y++)
		{
			adr = y*imgx+x;
		  /*
		     Norm 2 computation using 2x2 pixel window:
		       A B C
		       D E F
			   G H I
		     and
		       com1 = C-G,  com2 = I-A.
		     Then
		       gx = C+2F+I - (A+2D+G)=com1+com2+2(F-D)   horizontal difference
		       gy = G+2H+I - (A+2B+C)=-com1+com2+2(H-B)   vertical difference
		     com1 and com2 are just to avoid 2 additions.
		   */
		  com1 = img_in[adr-imgx+1] - img_in[adr+imgx-1];
		  com2 = img_in[adr+imgx+1] - img_in[adr-imgx-1];

		  gx = (com1+com2+2*(img_in[adr+1] - img_in[adr-1]))/(8.0*255); /* gradient x component */
		  gy = (-com1+com2+2*(img_in[adr+imgx] - img_in[adr-imgx]))/(8.0*255); /* gradient y component */
		  norm_square = gx*gx+gy*gy;
		  sum+=norm_square;

		  norm = sqrt( norm_square); /* gradient norm */

		  (mod)->data[adr] = norm; /* store gradient norm */
		   /* gradient angle computation */
	     (*angles)->data[adr] = atan2(gy,gx);
		}
	threshold = 2*sqrt(sum/(imgx*imgy));//自动阈值 (automatic threshold)
	//non maximum suppression
	for(x=1;x<imgx-1;x++)
		for(y=1;y<imgy-1;y++)
		{
			adr = y*imgx+x;
			value = (*angles)->data[adr];
			if((mod)->data[adr] < threshold )
			{
				(*angles)->data[adr] = NOTDEF;
				continue;
			}
			if( (value > -M_1_8_PI && value<=M_1_8_PI) || (value <= -M_7_8_PI ) || (value > M_7_8_PI))
			{
				if((mod)->data[adr] <= (mod)->data[adr+1] || (mod)->data[adr] <= (mod)->data[adr-1])
					(*angles)->data[adr] = NOTDEF;
			}
			else if( (value> M_1_8_PI && value<= M_3_8_PI) || (value> -M_7_8_PI && value<= -M_5_8_PI) )
			{
				if((mod)->data[adr] <= (mod)->data[adr-imgx-1] || (mod)->data[adr] <= (mod)->data[adr+imgx+1])
					(*angles)->data[adr] = NOTDEF;
			}
			else if((value> M_3_8_PI && value<= M_5_8_PI) || (value> -M_5_8_PI && value<= -M_3_8_PI))
			{
				if((mod)->data[adr] <= (mod)->data[adr-imgx] || (mod)->data[adr] <= (mod)->data[adr+imgx])
					(*angles)->data[adr] = NOTDEF;
			}
			else 
			{
				if((mod)->data[adr] <= (mod)->data[adr-imgx+1] || (mod)->data[adr] <= (mod)->data[adr+imgx-1])
					(*angles)->data[adr] = NOTDEF;
			}
		}
    //也标记到mod图上面 (Also marked on the mod map)
	//for(x=1;x<imgx-1;x++)
	//	for(y=1;y<imgy-1;y++)
	//	{
	//		if((*angles)->data[y*imgx+x] == NOTDEF)
	//			(mod)->data[y*imgx+x] = NOTDEF;
	//	}
		free_image_double(mod);
}

//canny
void calculateGradient3( double * img_in, unsigned int imgx, unsigned int imgy, image_double * angles)
{
	Mat1b edge;
	Mat1s DX,DY;
	Mat1b gray = Mat::zeros(imgy,imgx,CV_8UC1);
	unsigned int x,y,addr;
	(*angles) = new_image_double(imgx,imgy);
	//copy to gray image
	for ( y = 0; y<imgy; y++)
		for ( x = 0; x<imgx; x++)
		{
			addr = y*imgx+x;
			gray.data[addr] = (uchar)(img_in[addr]);
		}
	//canny
   Canny3(gray,edge,DX,DY,3,false);
   for ( y = 0; y<imgy; y++)
   {
	    short* _dx = DX.ptr<short>(y);
		short* _dy = DY.ptr<short>(y);
		uchar* _e = edge.ptr<uchar>(y);
		for ( x = 0; x<imgx; x++)
		{
			if(_e[x] > 0)//0 or 255
			{
				(*angles)->data[y*imgx+x]  = atan2((double)_dy[x],(double)_dx[x]);//calculate gradient 
			}
			else
				(*angles)->data[y*imgx+x] = NOTDEF;
		}
   }
   edge.release();
   DX.release();
   DY.release();
   gray.release();
}

#endif