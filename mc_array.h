/*
   MODULE MC_Array
   Implements dynamic array variables in one two and three Dimensions
*/

#ifndef INLUDE_MC_ARRAY
#define INLUDE_MC_ARRAY

#ifndef RANGECHECK
  #define INLINE inline
#else
  #undef INLINE
#endif

class MemoryHolder  {
  private:
    int *ncopy;
    void *base;
  public:
    MemoryHolder(void *buff);
    MemoryHolder(MemoryHolder *h);
    ~MemoryHolder(void);
};

template <class T> class Array1D  {
  private:
    int n, stride;
    T *data;
    MemoryHolder *base;
    void InRange(const int i);
  public:
    // Constructors
    Array1D<T>(const int n1);
    Array1D<T>(const int n1, T *buff, MemoryHolder *base, const int stride_in = 1);
    Array1D<T>(Array1D<T> &old);  	      // Not implemented!
    Array1D<T> operator = (Array1D<T> &old);  // Not implemented!
    ~Array1D(void);
    // Accessors
    int N(void)  { return n; }
    T & operator [] (const int i)  { return data[i]; }
    INLINE T & operator () (const int i);
};

template <class T> class Array2D  {
  private:
    int n[2], nn, stride[2];
    T *data;
    MemoryHolder *base;
    void InRange(const int i, const int j);
  public:
    // Constructors
    Array2D<T>(const int n1, const int n2);
    Array2D<T>(const int n1, const int n2, T *buff, MemoryHolder *base,
    	       const int stride_in[]);
    ~Array2D(void);
    Array2D<T>(Array2D<T> &old);  	      // Not implemented!
    Array2D<T> operator = (Array2D<T> &old);  // Not implemented!
    // Accessors
    int N(const int i)  { return n[i]; }
    T & operator [] (const int i)  { return data[i]; }
    INLINE T & operator () (const int i, const int j);
    // Util
    Array1D<T> *Slice(const int dim, const int i);
};

template <class T> class Array3D  {
  private:
    int n[3], nn, stride[3];
    T *data;
    MemoryHolder *base;
    void InRange(const int i, const int j, const int k);
  public:
    // Contructors
    Array3D<T>(const int n1, const int n2, const int n3);
    Array3D<T>(const int n1, const int n2, const int n3,
    	       T *buff, MemoryHolder *base, const int stride_in[]);
    ~Array3D(void);
    Array3D<T>(Array3D<T> &old);  	      // Not implemented!
    Array3D<T> operator = (Array3D<T> &old);  // Not implemented!
    // Accessors
    int N(const int i)  { return n[i]; }
    T & operator [] (const int i)  { return data[i]; }
    INLINE T & operator () (const int i, const int j, const int k);
    // Util
    Array1D<T> *Slice(const int dim, const int i, const int j);
    Array2D<T> *Slice(const int dim, const int i);
};

#ifndef RANGECHECK

template <class T>
T & Array1D<T>::operator () (const int i)
{
  return data[i*stride];
}

template <class T>
T & Array2D<T>::operator () (const int i, const int j)
{
  return data[i*stride[0]+j*stride[1]];
}

template <class T>
T & Array3D<T>::operator () (const int i, const int j, const int k)
{
  return data[i*stride[0]+j*stride[1]+k*stride[2]];
}
#endif

#endif
