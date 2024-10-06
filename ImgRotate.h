////////////////////////////////////////////////////////////////////////////////
//
// Fast and high quality true-color bitmap rotation function.
//
// Copyright (c) 2003 by Anton Treskunov
// Permission to use, copy, modify, distribute and sell this software for any 
//     purpose is hereby granted without fee, provided that the above copyright 
//     notice appear in all copies and that both that copyright notice and this 
//     permission notice appear in supporting documentation.
// The author makes no representations about the 
//     suitability of this software for any purpose. It is provided "as is" 
//     without express or implied warranty.
//
// The code is Pentium MMX optimized implementation of bitmap rotation by means 
//     of 3 shears.
// For the algorithm see: Alan Paeth, "A Fast Algorithm for General Raster 
//     Rotation," Graphics Interface '86, pp. 77-81, May 1986. For online
//     explanation refer to Tobin Fricke notes 
//     http://splorg.org/people/tobin/projects/israel/projects/paeth/rotation_by_shearing.html
// The software was built from "High quality image rotation (rotate by shear)", 
//     by Eran Yariv function on codeproject:
//     http://www.codeproject.com/bitmap/rotatebyshear.asp
// This implementation is about 6 times faster (tested on Athlon 1.4) and
//     more flexible.
// Speed is achived by 1) using byte weight factor instead of double; 
//     2) implementing internal shear loop in assembler using MMX instructions.
// Flexibility is achived by using design policies, inspired by 
//     "Modern C++ Design", by Andrei Alexandrescu, http://www.moderncppdesign.com
//
//	Class ImageRotate, parameterized by callback type, Allocator and Alignment
//     policies, contains 
//     	class Img	with members:
//			Pixel * pixel;		// struct { unsigned __int8 c[4]; }
//			unsigned width;		// in Pixels
//			unsigned height;	// in Pixels
//			unsigned width_pad;	// in bytes
//     and 
//      static Img AllocAndRotate (  
//			const Img &src,
//			Pixel clrBack,			// Background color
//			double      angle,      // Rotation angle
//			ProgressAnbAbortCallBack *cb = 0);
//  Callback type should support [bool operator (double)] and don't raise exceptions.
//     It could be function pointer (default) or functor class (e.g. Loki::Functor)
//	   Function gets double value of progress percentage in 0-100 range and 
//     returns false to abort calculation.
//  AllocatorPolicy should provide 2 methods w/o exceptions raising:
//	   void* allocate(size_t) and void free(void *t).
//     It could be "C" malloc/free (default) or custom memory manager.
//  AlignmentPolicy should provide 1 method w/o exception raising:
//     	unsigned aligned_width(unsigned width). By default, it implements
//      Windows BMP alignment.
//
//  Source and rotated images are kept in ImageRotate::Img classes.
//  Client is responsible to delete returned rotated image pixel array,
//     using image.Free();
//
////////////////////////////////////////////////////////////////////////////////

// Last update: March 20, 2003

#ifndef IMAGEROTATE_H
#define IMAGEROTATE_H

typedef  struct { unsigned __int8 c[4]; }   Pixel;

inline Pixel operator - (const Pixel a, Pixel b)
{
	Pixel d (a);
	d.c[0] -= b.c[0];
	d.c[1] -= b.c[1];
	d.c[2] -= b.c[2];
	d.c[3] -= b.c[3];
	return d;
}

inline Pixel operator * (const Pixel a, unsigned __int8 w)
{
	Pixel d;
	d.c[0] = (int)a.c[0] * w / 256;
	d.c[1] = (int)a.c[1] * w / 256;
	d.c[2] = (int)a.c[2] * w / 256;
	d.c[3] = (int)a.c[3] * w / 256;
	return d;
}

// AllocatePolicy

class AllocatePolicyStdNew
{
public:
	static void* allocate(size_t size) 
	{ return std::malloc(size); 	}
	static void free(void *t) { std::free (t); }
protected:
	~AllocatePolicyStdNew(){} // to prohibit destruction by client
};

#ifdef __GC		// Boehm Garbage Collector
class AllocatePolicyStdGC
{
public:
	static void* allocate(size_t size) 
	{ return GC_malloc_atomic(size); 	}
	static void free(void *t) {}
protected:
	~AllocatePolicyStdGC(){} // to prohibit destruction by client
};
#endif

// AlignmentPolicy

class AlignmentPolicyBmp
{
public:
	static unsigned aligned_width(unsigned width) { return (sizeof(Pixel)*width + 3) & (~3); }
protected:
	~AlignmentPolicyBmp(){} // to prohibit destruction by client
};

class AlignmentPolicyUnaligned
{
public:
	static unsigned aligned_width(unsigned width) { return sizeof(Pixel)*width; }
protected:
	~AlignmentPolicyUnaligned(){} // to prohibit destruction by client
};

// Callback

typedef bool (*CallbackFn)(double percent_complete);
#ifdef USE_LOKI
typedef Loki::Functor<bool, TYPELIST_1(double)> CallbackFtr;
#endif

template<
	typename ProgressAnbAbortCallBack = CallbackFn,
	class AllocatorPolicy = AllocatePolicyStdNew,
	class AlignmentPolicy = AlignmentPolicyBmp
>
class ImageRotate : 
	public AllocatorPolicy
{
public:
	class Img :
		public AllocatorPolicy,
		public AlignmentPolicy
	{ 
	public:
		Img () :
			width(0), height(0), width_pad(0), pixel(0)
		{}
		Img (const Img &s) :
			width(s.width), height(s.height), width_pad(s.width_pad), pixel(s.pixel)
		{}
		Img (int w, int h)
		{
			Allocate (w, h);
		}
		void Allocate (int w, int h)
		{
			width = w;
			height = h;
			width_pad = aligned_width(width);
			pixel = (Pixel *)allocate (width_pad * height);
		}
		void Free ()
		{
			free (pixel);
			pixel = 0;
		}

		Pixel * pixel; 
		unsigned width;
		unsigned height; 
		unsigned width_pad;	// in bytes

		Pixel &RGBVAL(int x, int y)
		{
			return *((Pixel *)((char *)pixel + width_pad * y) + x);
		}

		const Pixel &RGBVAL(int x, int y) const
		{
			return *((Pixel *)((char *)pixel + width_pad * y) + x);
		}

		void PrevLine (Pixel *& p) const
		{
			p = (Pixel *) ((char *)p - width_pad);
		}

		void NextLine (Pixel *& p) const
		{
			p = (Pixel *) ((char *)p + width_pad);
		}
	};

	static
	Img AllocAndRotate (  
		const Img &src,
		Pixel clrBack,			// Background color
		double      angle,      // Rotation angle
		ProgressAnbAbortCallBack *cb = 0
		)
	{
		if (!src.pixel)
			return src;

		Img mid_image (src); 

		// Bring angle to range of (-INF .. 360.0)
		while (angle >= 360.0)
			angle -= 360.0;

		// Bring angle to range of [0.0 .. 360.0) 
		while (angle < -0.0)
			angle += 360.0;

		if (angle > 45.0 && angle <= 135.0) 
		{
			// Angle in (45.0 .. 135.0] 
			// Rotate image by 90 degrees into temporary image,
			// so it requires only an extra rotation angle 
			// of -45.0 .. +45.0 to complete rotation.
			mid_image = Rotate90 (src, cb);
			angle -= 90.0;
		}
		else if (angle > 135.0 && angle <= 225.0) 
		{ 
			// Angle in (135.0 .. 225.0] 
			// Rotate image by 180 degrees into temporary image,
			// so it requires only an extra rotation angle 
			// of -45.0 .. +45.0 to complete rotation.
			mid_image = Rotate180 (src, cb);
			angle -= 180.0;
		}
		else if (angle > 225.0 && angle <= 315.0) 
		{ 
			// Angle in (225.0 .. 315.0] 
			// Rotate image by 270 degrees into temporary image,
			// so it requires only an extra rotation angle 
			// of -45.0 .. +45.0 to complete rotation.
			mid_image = Rotate270 (src, cb);
			angle -= 270.0;
		}

		// check for abort
		if (!mid_image.pixel)
			return mid_image;

		// If we got here, angle is in (-45.0 .. +45.0]

		Img dst (Rotate45 (mid_image, clrBack, angle, src.pixel != mid_image.pixel, cb));

		if (src.pixel != mid_image.pixel)
		{
			// Middle image was required, free it now.
			mid_image.Free();
		}

		return dst;
	}
	private:

	static
	Img Rotate90  (const Img &src, ProgressAnbAbortCallBack *cb = 0)
	{
		Img dst (src.height, src.width);
		if (!dst.pixel)
			return dst;

		for (unsigned uY = 0; uY < src.height; uY++)
		{
			const Pixel *psrc = &src.RGBVAL(0, uY);
			Pixel *pdst = &dst.RGBVAL (uY, dst.height-1);
			for (unsigned uX = 0; uX < src.width; uX++)
			{
				*pdst = *psrc;
				psrc++;
				dst.PrevLine(pdst);
			}
			// Report progress
			if (cb && !((*cb) (50.0 * uY / src.height)))
			{
				dst.Free();
				break;
			}
		}
		return dst;
	} 

	static
	Img Rotate180  (const Img &src, ProgressAnbAbortCallBack *cb = 0)
	{
		Img dst (src.width, src.height);
		if (!dst.pixel)
			return dst;

		for (unsigned uY = 0; uY < src.height; uY++)
		{
			const Pixel *psrc = &src.RGBVAL(0, uY);
			Pixel *pdst = &dst.RGBVAL (dst.width -1, dst.height - uY - 1);
			for (unsigned uX = 0; uX < src.width; uX++)
			{
				*pdst = *psrc;
				psrc++;
				pdst--;
			}
			// Report progress
			if (cb && !((*cb) (50.0 * uY / src.height)))
			{
				dst.Free();
				break;
			}
		}
		return dst;
	}

	static
	Img Rotate270  (const Img &src, ProgressAnbAbortCallBack *cb = 0)
	{
		Img dst (src.height, src.width);
		if (!dst.pixel)
			return dst;

		for (unsigned uY = 0; uY < src.height; uY++)
		{
			const Pixel *psrc = &src.RGBVAL(0, uY);
			Pixel *pdst = &dst.RGBVAL (dst.width - uY - 1, 0);
			for (unsigned uX = 0; uX < src.width; uX++)
			{
				*pdst = *psrc;
				psrc++;
				dst.NextLine(pdst);
			}
			// Report progress
			if (cb && !((*cb) (50.0 * uY / src.height)))
			{
				dst.Free();
				break;
			}
		}
		return dst;
	} 
	
	static
	Img Rotate45  (
		const Img &src,		// Img to rotate (+dimensions)
		Pixel clrBack,			// Background color
		double dAngle,			// Degree of rotation
		bool bMidImage,			// Was middle image used (for correct progress report)
		ProgressAnbAbortCallBack *cb = 0
		)
	{
		const double ROTATE_PI = 3.1415926535897932384626433832795;
		double dRadAngle = dAngle * ROTATE_PI / double(180); // Angle in radians
		double dSinE = sin (dRadAngle);
		double dTan = tan (dRadAngle / 2.0);

		// Calc first shear (horizontal) destination image dimensions 
		Img dst1 (src.width + int(double(src.height) * fabs(dTan) + 0.5), src.height);
		if (!dst1.pixel)
			return dst1;


		// Perform 1st shear (horizontal) 
		for (unsigned u = 0; u < dst1.height; u++) 
		{  
			double dShear = ((dTan >= 0.0 ? (int)u : (int)u - (int)dst1.height) + 0.5) * dTan ;
			int iShear = (int)floor (dShear);
			HorizSkew (src, dst1, u, iShear, unsigned __int8(255 * (dShear - double(iShear)) + 1), clrBack);
			// Report progress
			if (cb && !((*cb) (
				bMidImage ? 50.0 + (50.0 / 3) * u / dst1.height :
			                (100.0 / 3) * u / dst1.height
				)))
			{
				dst1.Free();
				return dst1;
			}
		}

		// Perform 2nd shear  (vertical)
		Img dst2 (dst1.width, int(double (src.width) * fabs (dSinE) + double (src.height) * cos (dRadAngle) + 0.5) + 1);
		if (!dst2.pixel)
		{
			dst1.Free();
			return dst2;
		}
		double dOffset     // Variable skew offset
			= (dSinE > 0.0 ? (int)src.width - 1 : (int)dst2.width - (int)src.width) * dSinE; 
		for (u = 0; u < dst2.width; u++, dOffset -= dSinE) 
		{
			int iShear = int (floor (dOffset));
			VertSkew (dst1, dst2, u, iShear, unsigned __int8(255 * (dOffset - double(iShear)) + 1), clrBack);
			// Report progress
			if (cb && !((*cb) (
				bMidImage ? 66.0 + (50.0 / 3) * u / dst2.height :
			                33.0 + (100.0 / 3) * u / dst2.height
				)))
			{
				dst1.Free();
				dst2.Free();
				return dst2;
			}
		}
		// Free result of 1st shear    
		dst1.Free();

		// Perform 3rd shear (horizontal) 
		Img dst3 (int (double(src.height) * fabs (dSinE) + double(src.width) * cos (dRadAngle)) + 1, dst2.height);
		if (!dst3.pixel)
		{
			dst2.Free();
			return dst3;
		}
		dOffset = dSinE >= 0.0 ? 
			(src.width - 1) * dSinE * -dTan 
			:
			dTan * (double(src.width - 1) * -dSinE + double(1 - (int)dst3.height));
		for (u = 0; u < dst3.height; u++, dOffset += dTan)
		{
			int iShear = int (floor(dOffset));
			HorizSkew (dst2, dst3, u, iShear, unsigned __int8(255 * (dOffset - double (iShear)) + 1), clrBack);
			if (cb && !((*cb) (
				bMidImage ? 83.0 + (50.0 / 3) * u / dst3.height :
			                66.0 + (100.0 / 3) * u / dst3.height
				)))
			{
				dst2.Free();
				dst3.Free();
				return dst3;
			}
		}
		// Free result of 2nd shear    
		dst2.Free();

		return dst3;      
	} 

	static
	void HorizSkew (
		const Img &src,		// Img to skew (+ dimensions)
		Img &dst,				// Destination of skewed image (+ dimensions)
		unsigned uRow,			// Row index
		int iOffset,			// Skew offset 
		unsigned __int8 Weight,			// Relative weight of right pixel
		Pixel clrBack			// Background color
		)
	{
		// Fill gap left of skew with background
		Pixel *p = &dst.RGBVAL(0, uRow);
		for (int i = 0u; i < iOffset; i++)
			*p++ = clrBack;


		Pixel pxlOldLeft;
		if (MMXpresnt())
		{	
			int iw = (Weight<<24)|(Weight<<16)|(Weight<<8)|Weight;
			int iw1 = ~iw & 0xFFFFFFFF;
			short *sptr = (short *)&src.RGBVAL (0, uRow);
			short *dptr = (short *)&dst.RGBVAL (iOffset, uRow);
			int swidth = src.width;
			int dwidth = dst.width;

			_asm pxor mm0, mm0	// mm0 -zero
			_asm movd mm1, clrBack	
			_asm PUNPCKLBW mm1, mm0 // mm1 - back
			_asm movq mm5, mm1      // mm5 - pxlOldLeft
			_asm mov edx, iw
			_asm movd mm2, edx	
			_asm PUNPCKLBW mm2, mm0 // mm2 - weight
			_asm mov edx, iw1
			_asm movd mm3, edx	
			_asm PUNPCKLBW mm3, mm0 // mm3 - 1-weight
			_asm PMULLW mm3, mm1	// mm3 - back*(1-weight)
			_asm mov esi, dword ptr[sptr]	// src
			_asm mov edi, dword ptr[dptr]	// dst

			_asm xor ecx, ecx		// ecx - loop counter
			_asm mov ebx, iOffset	// ebx - i+iOffset
	loop:   _asm cmp ecx, swidth
			_asm jge endloop

			_asm MOVD mm6, word ptr[esi]
			_asm PUNPCKLBW mm6, mm0	// mm6 - pxlSrc
			_asm movq mm4, mm6   	// mm4 - pixel
			_asm PMULLW mm4, mm2	// mm4 - pixel*weight
			_asm PADDW mm4, mm3		// mm4 - pixel*weight + back*(1-weight)
			_asm PSRLW mm4, 8		// mm4 - (pixel*weight + back*(1-weight)) / 256 - pxlLeft

			_asm cmp ebx, 0			// if ((i + iOffset >= 0) && (i + iOffset < dst.width)) 
			_asm jl noupdate
			_asm cmp ebx, dwidth
			_asm jge noupdate

			// dst.RGBVAL (i + iOffset, uRow) = pxlSrc - (pxlLeft - pxlOldLeft);
			_asm PSUBW mm6, mm4
			_asm PADDW mm6, mm5
			_asm PACKUSWB mm6, mm0	// mm6 - packed (pxlSrc - (pxlLeft - pxlOldLeft))
			_asm MOVD word ptr[edi], mm6

	noupdate:						// end of if

			_asm MOVQ mm5, mm4		// pxlOldLeft = pxlLeft

			_asm inc ebx			// i+iOffset
			_asm add esi, 4
			_asm add edi, 4

			_asm inc ecx			// end of for
			goto loop;
	endloop:
			_asm PACKUSWB mm5, mm0
			_asm MOVD pxlOldLeft, mm5			// pxlOldLeft
			_asm EMMS							// restore FPU status

			// update i, iOffset, pxlOldLeft
		}
		else
		{
			pxlOldLeft = clrBack;

			for (unsigned i = 0; i < src.width; i++) 
			{
				// Loop through row pixels
				Pixel pxlSrc = src.RGBVAL (i, uRow);
				// Calculate weights
				Pixel pxlLeft = interp (pxlSrc, clrBack, Weight);
				// Check boundries 
				if ((i + iOffset >= 0) && (i + iOffset < dst.width))
				{
					// Update left over on source
					dst.RGBVAL (i + iOffset, uRow) = pxlSrc - (pxlLeft - pxlOldLeft);
				}
				// Save leftover for next pixel in scan
				pxlOldLeft = pxlLeft;
			}
		}
		// Go to rightmost point of skew
		i = src.width + iOffset;  
		// If still in image bounds, put leftovers there
		if (i < dst.width)
			dst.RGBVAL (i, uRow) = pxlOldLeft;
		p = &dst.RGBVAL(i, uRow);
		while (++i < dst.width)
		{   // Clear to the right of the skewed line with background
			*++p = clrBack;
		}
	}

	// Skews a column vertically (with filtered weights)
	// Limited to 45 degree skewing only. Filters two adjacent pixels.
	static
	void VertSkew (
		const Img &src,		// Img to skew (+ dimensions)
		Img &dst,				// Destination of skewed image (+ dimensions)
		unsigned uCol,				// Column index
		int iOffset,			// Skew offset 
		unsigned __int8 Weight,			// Relative weight of right pixel
		Pixel clrBack			// Background color
		)
	{
		Pixel *p = &dst.RGBVAL(uCol, 0);
		for (int i = 0; i < iOffset; i++)
		{
			// Fill gap above skew with background
			*p = clrBack;
			dst.NextLine (p);
		}

		Pixel pxlOldLeft;
		if (MMXpresnt())
		{	
			int iw = (Weight<<24)|(Weight<<16)|(Weight<<8)|Weight;
			int iw1 = ~iw & 0xFFFFFFFF;
			short *sptr = (short *)&src.RGBVAL (uCol, 0);
			short *dptr = (short *)&dst.RGBVAL (uCol, iOffset);
			int sheight = src.height;
			int dheight = dst.height;
			int swidthpad = src.width_pad;
			int dwidthpad = dst.width_pad;

			_asm pxor mm0, mm0	// mm0 -zero
			_asm movd mm1, clrBack	
			_asm PUNPCKLBW mm1, mm0 // mm1 - back
			_asm movq mm5, mm1      // mm5 - pxlOldLeft
			_asm mov edx, iw
			_asm movd mm2, edx	
			_asm PUNPCKLBW mm2, mm0 // mm2 - weight
			_asm mov edx, iw1
			_asm movd mm3, edx	
			_asm PUNPCKLBW mm3, mm0 // mm3 - 1-weight
			_asm PMULLW mm3, mm1	// mm3 - back*(1-weight)
			_asm mov esi, dword ptr[sptr]	// src
			_asm mov edi, dword ptr[dptr]	// dst

			_asm xor ecx, ecx		// ecx - loop counter
			_asm mov ebx, iOffset	// ebx - i+iOffset
	loop:   _asm cmp ecx, sheight
			_asm jge endloop

			//_asm mov ax, word ptr[esi]
			//_asm MOVD mm6, eax
			_asm MOVD mm6, word ptr[esi]
			_asm PUNPCKLBW mm6, mm0	// mm6 - pxlSrc
			_asm movq mm4, mm6   	// mm4 - pixel
			_asm PMULLW mm4, mm2	// mm4 - pixel*weight
			_asm PADDW mm4, mm3		// mm4 - pixel*weight + back*(1-weight)
			_asm PSRLW mm4, 8		// mm4 - (pixel*weight + back*(1-weight)) / 256 - pxlLeft

			_asm cmp ebx, 0			// if ((i + iOffset >= 0) && (i + iOffset < dst.height)) 
			_asm jl noupdate
			_asm cmp ebx, dheight
			_asm jge noupdate

			// dst.RGBVAL (i + iOffset, uRow) = pxlSrc - (pxlLeft - pxlOldLeft);
			_asm PSUBW mm6, mm4
			_asm PADDW mm6, mm5
			_asm PACKUSWB mm6, mm0	// mm6 - packed (pxlSrc - (pxlLeft - pxlOldLeft))
			_asm MOVD word ptr[edi], mm6

	noupdate:						// end of if

			_asm MOVQ mm5, mm4		// pxlOldLeft = pxlLeft

			_asm inc ebx			// i+iOffset
			_asm add esi, swidthpad
			_asm add edi, dwidthpad

			_asm inc ecx			// end of for
			goto loop;
	endloop:
			_asm PACKUSWB mm5, mm0
			_asm MOVD pxlOldLeft, mm5			// pxlOldLeft
			_asm EMMS							// restore FPU status

			// update i, iOffset, pxlOldLeft
		}
		else
		{
			pxlOldLeft = clrBack;

			for (i = 0; i < src.height; i++) 
			{
				// Loop through column pixels
				Pixel pxlSrc = src.RGBVAL (uCol, i);
				// Calculate weights
				Pixel pxlLeft = interp (pxlSrc, clrBack, Weight);
				// Update left over on source
				// Check boundries
				int iYPos = i + iOffset;
				if ((iYPos >= 0) && (iYPos < dst.height))
					dst.RGBVAL (uCol, iYPos) = pxlSrc - (pxlLeft - pxlOldLeft);
				// Save leftover for next pixel in scan
				pxlOldLeft = pxlLeft;
			}
		}
		// Go to bottom point of skew
		i = src.height + iOffset;  
		// If still in image bounds, put leftovers there
		if (i < dst.height)
			dst.RGBVAL (uCol, i) = pxlOldLeft;
		// Clear below skewed line with background
		p = &dst.RGBVAL(uCol, i);
		while (++i < dst.height)
		{
			dst.NextLine (p);
			*p = clrBack;
		}
	}
	
	static
	bool MMXpresnt()
	{
		return false;
		__asm mov EAX,1              ;request CPU feature flags
		__asm cpuid                  ;0Fh, 0A2h CPUID instruction
		__asm test EDX,800000h       ;test bit 23 to see if MMX available
		__asm jnz yes	             ;yes
		return false;
	yes:
		return true;
	}

	static
	unsigned __int8 interp(unsigned __int8 a, unsigned __int8 b, unsigned __int8 w)
	{ return b + (a-b)*w/256; }

	static
	Pixel interp(Pixel a, Pixel b, unsigned __int8 w)
	{
		Pixel dst;
		dst.c[0] = interp (a.c[0], b.c[0], w);
		dst.c[1] = interp (a.c[1], b.c[1], w);
		dst.c[2] = interp (a.c[2], b.c[2], w);
		dst.c[3] = interp (a.c[3], b.c[3], w);
		return dst;
	}

};

#endif
