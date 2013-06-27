/*
   MODULE MC_Array
   Implements dynamic array variables in one two and three Dimensions
   : template function definitions
*/

template <class T>
Array1D<T>::Array1D(const int n1)
{
  data = new T[n1];
  if (!data)  {
    cerr << "Error allocating " << n1
         << " elements in Array1D<T>::Array1D...aborting" << endl;
    abort();
  }
  base = new MemoryHolder(data);
  n = n1;
  stride = 1;
}

template <class T>
Array1D<T>::Array1D(const int n1, T *buff, MemoryHolder *base_in, const int stride_in = 1)
{
  data = buff; base = new MemoryHolder(base_in);
  n = n1;
  stride = stride_in;
}

template <class T>
Array1D<T>::~Array1D(void)
{
  delete base;
}

template <class T>
void Array1D<T>::InRange(const int i)
{
  if (i < 0 || i >= n)  {
    cerr << "Implementation ERROR: Index out of range (0.." << n-1
    	 << ") is " << i << endl;
    abort();
  }
}

template <class T>
Array2D<T>::Array2D(const int n1, const int n2)
{
  nn = n1 * n2;
  data = new T[nn];
  if (!data)  {
    cerr << "Error allocating " << nn
         << " elements in Array2D<T>::Array2D...aborting" << endl;
    abort();
  }
  n[0] = n1; n[1] = n2;
  stride[0] = 1; stride[1] = n1;
}

template <class T>
Array2D<T>::Array2D(const int n1, const int n2, T *buff,
		    MemoryHolder *base_in, const int stride_in[])
{
  nn = n1 * n2;
  data = buff;
  base = new MemoryHolder(base_in);
  n[0] = n1; n[1] = n2;
  stride[0] = stride_in[0]; stride[1] = stride_in[1];
}

template <class T>
Array2D<T>::~Array2D(void)
{
  delete base;
}

template <class T>
Array1D<T> *Array2D<T>::Slice(const int dim, const int i)
{
#ifdef RANGECHECK
  if (dim < 0 || dim > 1)  {
    cerr << "Implementation ERROR: Wrong Dimension number: " << dim << "!" << endl;
    exit (3);
  }
  if (dim)  InRange(0, i);
  else	    InRange(i, 0);
#endif
  return new Array1D<T>(n[1-dim], data+i*stride[dim], base, stride[1-dim]);
}

template <class T>
void Array2D<T>::InRange(const int i, const int j)
{
  if (i < 0 || i >= n[0] || j < 0 || j >= n[1])  {
    cerr << "Implementation ERROR: Index out of range i(0.." << n[0]-1
    	 << ") is " << i
    	 << "; j(0.." << n[1]-1 << ") is " << j << endl;
    abort();
  }
}

template <class T>
Array3D<T>::Array3D(const int n1, const int n2, const int n3)
{
  nn = n1 * n2 * n3;
  data = new T[nn];
  if (!data)  {
    cerr << "Error allocating " << nn
         << " elements in Array2D<T>::Array2D...aborting" << endl;
    abort();
  }
  n[0] = n1; n[1] = n2; n[2] = n3;
  stride[0] = 1; stride[1] = n1; stride[2] = n1*n2;
}

template <class T>
Array3D<T>::Array3D(const int n1, const int n2, const int n3, T *buff,
		    MemoryHolder *base_in, const int stride_in[])
{
  nn = n1 * n2 * n3;
  data = buff;
  base = new MemoryHolder(base_in);
  n[0] = n1; n[1] = n2; n[2] = n3;
  stride[0] = stride_in[0]; stride[1] = stride_in[1]; stride[2] = stride_in[2];
}

template <class T>
Array3D<T>::~Array3D(void)
{
  delete base;
}

template <class T>
Array1D<T> * Array3D<T>::Slice(const int dim, const int i, const int j)
{
#ifdef RANGECHECK
  switch (dim)  {
    case 0 : InRange(0, i, j); break;
    case 1 : InRange(i, 0, j); break;
    case 2 : InRange(i, j, 0); break;
  }
#endif
  switch (dim)  {
    case 0 : return new Array1D<T>(n[0], data + i*stride[1]+j*stride[2], base, stride[0]);
    case 1 : return new Array1D<T>(n[1], data + i*stride[0]+j*stride[2], base, stride[1]);
    case 2 : return new Array1D<T>(n[2], data + i*stride[0]+j*stride[1], base, stride[2]);
    default : cerr << "Implementation Error: Dimension in Array3D<T>::Slice is " << dim << endl;
  }
}

template <class T>
Array2D<T> * Array3D<T>::Slice(const int dim, const int i)
{
  int newstride[2];
#ifdef RANGECHECK
  switch (dim)  {
    case 0 : InRange(i, 0, 0); break;
    case 1 : InRange(0, i, 0); break;
    case 2 : InRange(0, 0, i); break;
  }
#endif
  switch (dim)  {
    case 0 : return new Array2D<T>(n[1], n[2], data + i*stride[0], base, stride+1);
    case 1 : newstride[0] = stride[0]; newstride[1] = stride[2];
    	     return new Array2D<T>(n[0], n[2], data + i*stride[1], base, newstride);
    case 2 : return new Array2D<T>(n[0], n[1], data + i*stride[2], base, stride);
    default : cerr << "Implementation Error: Dimension in Array3D<T>::Slice is " << dim << endl;
  }
}

template <class T>
void Array3D<T>::InRange(const int i, const int j, const int k)
{
  if (i < 0 || i >= n[0] || j < 0 || j >= n[1] || k < 0 || k >= n[2])  {
    cerr << "Implementation ERROR: Index out of range i(0.." << n[0]-1
    	 << ") is " << i
    	 << "; j(0.." << n[1]-1 << ") is " << j
    	 << "; k(0.." << n[2]-1 << ") is " << k << endl;
    abort();
  }
}

#ifdef RANGECHECK

template <class T>
T & Array1D<T>::operator () (const int i)
{
  InRange(i);
  return data[i*stride];
}

template <class T>
T & Array2D<T>::operator () (const int i, const int j)
{
  InRange(i, j);
  return data[i*stride[0]+j*stride[1]];
}

template <class T>
T & Array3D<T>::operator [] (const int i)
{
  return data[i];
}

template <class T>
T & Array3D<T>::operator () (const int i, const int j, const int k)
{
  InRange(i, j, k);
  return data[i*stride[0]+j*stride[1]+k*stride[2]];
}

#endif

