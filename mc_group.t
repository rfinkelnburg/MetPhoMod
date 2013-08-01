#include <cstddef>
#include <iostream>
#include <mc_group.hh>

template <class T> class Holder {
public :
  T *v;
  int id;
  Holder<T> *next;
};

template <class T>
Group<T>::Group(void)
{
  first = NULL;
  maxid = 0;
}

template <class T>
Group<T>::~Group(void)
{
  Holder<T> *vh, *h;
  for (vh = first; vh; ) {
    h = vh;
    vh = vh->next;
    delete h;
  }
}

template <class T>
int Group<T>::Register(T *v)
{
  Holder<T> *vh, *p;
  if (first) {
    for (p = first; p->next; p = p->next) {
      if (p->v == v)  // Item already registered, return its ID
	return p->id;
    }
    p->next = vh = new Holder<T>;
  }
  else {
    first = vh = new Holder<T>;
  }
  vh->v = v;
  vh->id = ++maxid;
  vh->next = NULL;
  return vh->id;
}

template <class T>
void Group<T>::Unregister(const T *v)
{
  Holder<T> *vh, *h;
  if (first->v == v) {
    h = first;
    first = first->next;
    delete h;
  }
  else {
    for (vh = first; vh->next && vh->next->v != v; vh = vh->next);
    if (vh->next) {
      h = vh->next;
      vh->next = vh->next->next;
      delete h;
    }
    else {
      std::cerr << "Implementation error: Variable can not be unregistered!\n";
      exit (3);
    }
  }
}

template <class T>
short Group<T>::Contains(const T *v) const
{
  Holder<T> *vh;
  for (vh = first; vh && vh->v != v; vh = vh->next);
  return (!!vh);
}

template <class T>
T *Group<T>::Find(const char *name) const
{
  Holder<T> *p;
  for (p = first; p && strcmp(name, p->v->name); p = p->next);
  if (p)  return (p->v);
  else    return (NULL);
}

template <class T>
T *Group<T>::Find(const int id) const
{
  Holder<T> *p;
  for (p = first; p && p->id != id; p = p->next);
  if (p) return (p->v);
  else   return (NULL);
}

template <class T>
int Group<T>::N(void) const
{
  Holder<T> *p;
  int n = 0;
  for (p = first; p; p = p->next) n++;
  return n;
}

template <class T>
GroupEnumerator<T>::GroupEnumerator(const Group<T> &group)
{
  p = group.first;
}

template <class T>
void GroupEnumerator<T>::SetGroup(const Group<T> &group)
{
  p = group.first;
}

template <class T>
short GroupEnumerator<T>::operator ! (void)
{
  return !!p;
}

template <class T>
T *GroupEnumerator<T>::operator () (void)
{
  return p->v;
}

template <class T>
void GroupEnumerator<T>::operator ++ (void)
{
  p = p->next;
}

template <class T>
int GroupEnumerator<T>::ID(void) const
{
  return p->id;
}
