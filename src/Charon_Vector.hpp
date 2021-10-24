
///////////////////////////////////////////////////////////////////////////////
//
//  Charon_Vector.hpp
//
///////////////////////////////////////////////////////////////////////////////
#ifndef   __Charon_Vector_hpp__
#define   __Charon_Vector_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <vector>

namespace charon
{

  /**
   *  \brief A Charon version of a `std::vector` with output streaming.
   *
   *  Normally when you `cout` a `ParameterList`, any `vector`s in the list get
   *  printed as memory addresses.  This class overload the output streaming
   *  operator such that you can actually see what's in those `vector`s.
   *
   *  \note This is simply useful for debugging.
   */
  template<typename T>
  class Vector
    :
    public std::vector<T>
  {
    public:

      /**
       *  \brief Default Constructor
       *
       *  Constructs an empty container with no elements.
       */
      Vector()
      {
      } // end of Default Constructor

      /**
       *  \brief Fill Constructor
       *
       *  Constructs a container with `n` elements.  Each element is a copy of
       *  a default-constructed object of type `T`.
       *
       *  \param[in] n The size of the container.
       */
      Vector(
        size_t n)
      {
        this->resize(n);
      } // end of Fill Constructor

      /**
       *  \brief Destructor.
       *
       *  This destroys all container elements and deallocates all the storage
       *  capacity allocated by the vector using its allocator.
       */
      ~Vector()
      {
        this->clear();
      } // end of Destructor

      /**
       *  \brief Output Streaming Operator
       *
       *  Output the contents of the `vector` to the given stream.
       *
       *  \param[in,out] os  The output stream to which you'd like to print.
       *  \param[in]     obj The vector you'd like to print to the stream.
       *
       *  \returns The modified output stream.
       */
      friend std::ostream& operator<<(
        std::ostream&  os,
        std::vector<T> obj)
      {
        if (obj.size() != 0)
          os << obj[0];
        return os;
      } // end of operator<<()

  }; // end of class Vector

} // end of namespace charon

#endif // __Charon_Vector_hpp__
