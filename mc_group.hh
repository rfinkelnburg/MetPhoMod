/*
   MODULE mc_group
   A group management template
*/

#ifndef INCLUDE_MC_GROUP
#define INCLUDE_MC_GROUP

class VariableManager;

template <class T> class Holder;
template <class T> class GroupEnumerator;

template <class T> class Group {
  Holder<T> *first;
  int maxid;
  friend class GroupEnumerator<T>;
public :
  // Creators 
  Group(void);
  ~Group(void);

  // Modifiers
  int Register(T *v);
  void Unregister(const T *v);

  // Accessors
  short Contains(const T *v) const;
  T *Find(const char *name) const;
  T *Find(const int id) const;
  void AddToManager(VariableManager &m) const;
  int N(void) const;
};

template <class T> class GroupEnumerator {
  Holder<T> *p;
public :
  // Creator
  GroupEnumerator(const Group<T> &group);
  
  // Modifiers
  void SetGroup(const Group<T> &group);
  short operator ! (void);  // Is there an active element?
  T * operator () (void); // The active element
  void operator ++ (void); // Switch to next element

  // Accessors
  int ID(void) const;  // Return the ID of the active element
};

#endif
