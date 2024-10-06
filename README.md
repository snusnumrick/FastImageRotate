# FastImageRotate

This is my old 2003 project, I put on GitHub now in hope it would be usefull to somebody.

## Sources and Motivation
This popular algorithm was proposed by Alan Paeth (*Alan Paeth, "A Fast Algorithm for General Raster Rotation," Graphics Interface '86, pp. 77-81, May 1986.*) For online explanation see [Tobin Fricke notes](#). The software implementation I came across on [codeproject.com](#) belongs to **Eran Yariv**. 

This implementation:
* while providing good quality is rather slow ("*Code is optimized for image quality, not speed.*")
* while claims to support "*generic bitmap structures*", always uses **COLORREF** as underlying pixel representation, which seriously limits flexibility.

## Benefits

My implementation is about 6 times faster (tested on Athlon 1.4) and more flexible.

Speed is achieved by:
* using byte weight factor instead of double for linear interpolation;
* implementing internal shear loop in assembler **using MMX instructions**;
* eliminating virtual functions for individual pixel access and using **static polymorphism** instead.

Flexibility is achieved by using design policies, inspired by **Andrei Alexandrescu**.


## Interface

Class ImageRotate is parameterized by callback type, Allocator and Alignment policies.

```cpp
template<
    typename ProgressAnbAbortCallBack = CallbackFn,
    class AllocatorPolicy = AllocatePolicyStdNew,
    class AlignmentPolicy = AlignmentPolicyBmp
>
class ImageRotate;
```

Source and destination images are described by **class ImageRotate::Img**

```cpp
Pixel *pixels
struct { unsigned __int8 c[4]; }
unsigned width in Pixels
unsigned height in Pixels
unsigned width_pad in bytes
```

Public interface of **ImageRotate** consists from single function

```cpp
static Img AllocAndRotate (  
    const Img &src,
    Pixel clrBack,     // Background color
    double      angle, // Rotation angle
    ProgressAnbAbortCallBack *cb = 0);
```

Callback type should support

```cpp
bool operator() (double) 
```

and may not raise exceptions. It could be function pointer (default) or functor class (e.g. **Loki::Functor**). Function gets double value of progress percentage in 0-100 range and returns *false* to abort calculation or *true* otherwise.

AllocatorPolicy should provide 2 methods w/o exceptions raising:

```cpp
void* allocate(size_t);
void free(void *t);
```

It could be "C" *malloc*/*free* (default) or custom memory manager.

AlignmentPolicy should provide 1 method w/o exception raising:

```cpp
unsigned aligned_width(unsigned width);
```

By default, it implements *Windows BMP* 4 byte alignment.

Source and rotated images are kept in ImageRotate::Img classes. Client is responsible to delete returned rotated image pixel array, using *image.Free()*;

## References

1. Alan Paeth, "A Fast Algorithm for General Raster Rotation," Graphics Interface '86, pp. 77-81, May 1986.
2. Tobin Fricke, [Rotation by Shearing](#)
3. Eran Yariv, [High quality image rotation (rotate by shear)](#)
4. Intel Corporation, [Using MMX™ Instructions to Implement Alpha Blending](#)
5. Todd Veldhuizen, [Techniques for Scientific C++](#)
6. Andrei Alexandrescu, [Modern C++ Design](#)
7. SourceFORGE.net, [Project: Loki](#)

Copyright © 2003 **Anton Treskunov**

Permission to copy, use, modify, sell and distribute this document is granted provided this copyright notice appears in all copies. This document is provided "as is" without express or implied warranty, and with no claim as to its suitability for any purpose.
