#ifndef MULI_PROMOTE_HXX
#define MULI_PROMOTE_HXX

namespace muli {

template <class T1, class T2 = T1>
using Promote = decltype(T1()+T2());

} // namespace muli

#endif // MULI_PROMOTE_HXX