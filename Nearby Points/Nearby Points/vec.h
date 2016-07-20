#ifndef DF_VEC
#define DF_VEC

#if defined(__SSE2__) || defined(_M_IX86_FP)
#include <emmintrin.h>
#endif

#include <cmath>
#include <iostream>

#pragma push_macro("DI")
#pragma push_macro("DEF_SQ_MAT_BLOCK_COMMON")
#pragma push_macro("DEF_MAT_SINGLE_OP")
#pragma push_macro("DEF_SIMPLE_MAT_MAT_OP")
#pragma push_macro("DEF_MAT_SINGLE_IDENTITY")
#pragma push_macro("LG")
#pragma push_macro("DEF_MAT_BLOCK_COMMON")
#pragma push_macro("DEF_MAT_SINGLE")
#pragma push_macro("DEF_MI")
#pragma push_macro("DEF_VEC")
#pragma push_macro("DEF_VEC_SINGLE")
#pragma push_macro("F")
#pragma push_macro("VI")
#pragma push_macro("DEF_MAT_BLOCK_E_WISE_OP_EQ_TILDA")
#pragma push_macro("DEF_VEC_COMMON")
#pragma push_macro("DEF_VEC_BLOCK")
#pragma push_macro("DEF_NAMESPACE")
#pragma push_macro("U")
#pragma push_macro("IL")
#pragma push_macro("DEF_SSE_INV")
#pragma push_macro("SX")
#pragma push_macro("DEF_MAT_BLOCK_E_WISE_OP_EQ")
#pragma push_macro("SS")
#pragma push_macro("DEF_MAT_TILDA_OP")
#pragma push_macro("MI")
#pragma push_macro("ST")
#pragma push_macro("DEF_MAT_SCALAR_OP")
#pragma push_macro("DEF_IS_ARI")
#pragma push_macro("SI")
#pragma push_macro("SM")
#pragma push_macro("DEF_SQ_MAT_COMMON")
#pragma push_macro("DEF_MAT_BLOCK_E_WISE_OP_EQ_SCALAR")
#pragma push_macro("DEF_MAT_SINGLE_OP_SCALAR")

#undef DI
#undef DEF_SQ_MAT_BLOCK_COMMON
#undef DEF_MAT_SINGLE_OP
#undef DEF_SIMPLE_MAT_MAT_OP
#undef DEF_MAT_SINGLE_IDENTITY
#undef LG
#undef DEF_MAT_BLOCK_COMMON
#undef DEF_MAT_SINGLE
#undef DEF_MI
#undef DEF_VEC
#undef DEF_VEC_SINGLE
#undef F
#undef VI
#undef DEF_MAT_BLOCK_E_WISE_OP_EQ_TILDA
#undef DEF_VEC_COMMON
#undef DEF_VEC_BLOCK
#undef DEF_NAMESPACE
#undef U
#undef IL
#undef DEF_SSE_INV
#undef SX
#undef DEF_MAT_BLOCK_E_WISE_OP_EQ
#undef SS
#undef DEF_MAT_TILDA_OP
#undef MI
#undef ST
#undef DEF_MAT_SCALAR_OP
#undef DEF_IS_ARI
#undef SI
#undef SM
#undef DEF_SQ_MAT_COMMON
#undef DEF_MAT_BLOCK_E_WISE_OP_EQ_SCALAR
#undef DEF_MAT_SINGLE_OP_SCALAR


#define DEF_NAMESPACE namespace DF {
DEF_NAMESPACE
#undef DEF_NAMESPACE

#define F(I, TILL) for(unsigned I = 0; I < TILL; ++I)
#define IL inline
#define U unsigned
#define VI(INDEX) this->operator[](INDEX)
#define MI(ROW, COL) this->operator()(ROW, COL)

#define DEF_MI \
IL S & operator()(U row, U col) \
{ return static_cast<Derived*>(this)->operator()(row, col); } \
IL S operator()(U row, U col) const \
{ return static_cast<const Derived*>(this)->operator()(row, col); }

#define LG(R, C) (R==1 ? C : C==1 ? R : 0)

// Shorthands
// R: NumRows, C: NumColumns, MP: MaxArgPos, P: Pos, N: Next
// MR: MaxRow, C: MaxColumn

template <class S, U Dims> struct Vec;
template <class S, U NumRows, U NumCols> struct Mat;

template <class S, class Derived, U MR, U MC, U R, U C> struct tMatBlock;
template <class S, U R, U C, bool IsIdentity> struct tMatSingle;
template <class S, class Derived, U R, U C> struct tVec2Base;
template <class S, class Derived, U R, U C> struct tVec3Base;
template <class S, class Derived, U R, U C> struct tVec4Base;
template <class S, class Derived, U R> struct tInvertibleSqMatBlock;
template <class S, class Derived, U R> struct tInvertibleSqMat;
template <class S, class Derived, U R> struct tSqMatBase;
template <class S, U R, U C, class Derived, bool IsRow, bool IsVec>
struct tMatBase;

namespace tMatUtil
{
    template <class T> struct IsAri { enum { value = false }; };
    
#define DEF_IS_ARI(T) template <> struct IsAri<T> { enum { value = true }; };
    
    DEF_IS_ARI(bool);
    DEF_IS_ARI(char);
    DEF_IS_ARI(unsigned char);
    DEF_IS_ARI(short);
    DEF_IS_ARI(unsigned short);
    DEF_IS_ARI(int);
    DEF_IS_ARI(unsigned int);
    DEF_IS_ARI(long);
    DEF_IS_ARI(unsigned long);
    DEF_IS_ARI(long long);
    DEF_IS_ARI(unsigned long long);
    DEF_IS_ARI(float);
    DEF_IS_ARI(double);
    DEF_IS_ARI(long double);
    
#undef DEF_IS_ARI
    
    template <typename T> class IsMatBase
    {
        struct Fallback { int DFMatBase; };
        struct Derived : T, Fallback {};
        
        template <typename Y, Y> struct Check;
        
        typedef char One[1];
        typedef char Two[2];
        
        template <typename Y>
        static One & func(Check<int Fallback::*, &Y::DFMatBase> *);
        
        template <typename Y> static Two & func(...);
    public:
        
        enum { value = sizeof(func<Derived>(0)) == 2 };
    };
    
    template <bool C, typename T = void> struct EnableIf { typedef T type; };
    template <typename T> struct EnableIf <false, T> { };
    
    template <U Dims> struct InitVec
    {
        template <class V, class T, bool IA> struct O
        { static IL void init(V &v, const T &t) { F(i, Dims) v[i] = t[i]; } };
        
        template <class V, class T> struct O <V, T, true>
        { static IL void init(V &v, const T &t) { F(i, Dims) v[i] = t; } };
        
        template <class V, class T> IL InitVec(V &v, const T &t)
        { O<V, T, IsAri<T>::value>::init(v, t); }
    };
    
    template <U R, U C> struct InitMat
    {
        // One SFINAE after another ... hehehe
        
        enum { D = R==1 ? C : C==1 ? R : 0 };
        
        template <class S, class T, class M> static IL
        void init3(S *_, const tMatBase<T, R, C, M, false, false> &m)
        {
            F(j, C) F(i, R) _[j*R+i] = m(i,j);
        }
        
        template <class S, class T> static IL
        typename tMatUtil::EnableIf<!tMatUtil::IsMatBase<T>::value>::type
        init2(S *_, const T &t)
        {
            F(j, C) F(i, R) _[j*R+i] = t(i,j);
        }
        
        template <class S, class T> static IL
        typename tMatUtil::EnableIf<tMatUtil::IsMatBase<T>::value>::type
        init2(S *_, const T &t)
        {
            init3(_, t);
        }
        
        template <class S, class T> static IL
        typename tMatUtil::EnableIf<!tMatUtil::IsAri<T>::value>::type
        init(S *_, const T &t)
        {
            init2(_, t);
        }
        
        template <class S, class T> static IL
        typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value>::type
        init(S *_, const T &t)
        {
            F(j, C) F(i, R) _[j*R+i] = (j == i ? t : 0);
        }
        
        template <class S, class T> IL InitMat(S *_, const T &t) { init(_, t); }
    };
    
    struct InitSingle
    {
        template <class S, class T, bool IA> struct O
        { static IL void init(S &s, const T &t) { s = t[0]; } };
        
        template <class S, class T> struct O <S, T, true>
        { static IL void init(S &s, const T &t) { s = t; } };
        
        template <class S, class T> IL InitSingle(S &s, const T &t)
        { O<S, T, IsAri<T>::value>::init(s, t); }
    };
};

// Hack to create new operators:
// *~ (element wise multiplication)
// /~ (element wise division)
template <class S, class Derived, U R, U C> struct tMatTilda
{
    const Derived &d;
    IL tMatTilda(const Derived &d): d(d) {}
    IL S operator()(U row, U col) const { return d(row, col); }
    
    friend std::ostream& operator << (std::ostream &os, const tMatTilda &m)
    {
        F(j, R) {
            os << (j == 0 ? '<' : ' ');
            F(i, C) os << ' ' << m(j,i);
            os << ' ' << (j == R - 1 ? '>' : '\n');
        }
        return os << std::flush;
    }
};

template <class S, U R, U C, class Derived, bool IsRow, bool IsVec>
struct tMatBase;

template <class M, U R, U C, U MP, U P>
struct tMatScalarAssigner
{
    M &m;
    IL tMatScalarAssigner(M &m): m(m) {}
    IL operator M &() { return m; }
    typedef tMatScalarAssigner<M, R, C, MP, P + 1> N;
    template <class T> IL N operator , (const T &t)
    { m(P/C, P%C) = t; return N(m); }
};

template <class M, U R, U C, U MP>
struct tMatScalarAssigner <M, R, C, MP, MP>
{
    M &m;
    IL tMatScalarAssigner(M &m): m(m) {}
    IL operator M &() { return m; }
    template <class T> IL M & operator , (const T &t)
    { m(R-1, C-1) = t; return m; }
};

template <class S, U R, U C, class Derived, bool IsRow, bool IsVec>
struct tMatBase
{
private:
    template <class M>
    IL M _convert() const { M m; F(j,C) F(i,R) m(i,j) = MI(i,j); return m; }
public:
    
    template <class V>
    IL operator V() const volatile
    { return const_cast<const tMatBase*>(this)->_convert<V>(); }
    
    DEF_MI;
    
    typedef int DFMatBase;
    
    typedef tMatScalarAssigner<tMatBase, R, C, R * C - (R * C > 0), 1>
    ScalarAssigner;
    
    template <class T>
    IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value,
    ScalarAssigner>::type operator << (const T &t)
    { this->operator()(0,0) = t; return ScalarAssigner(*this); }
    
    typedef tMatTilda<S, tMatBase, R, C> Tilda;
    IL Tilda operator ~() const { return Tilda(*this); }
    
    friend std::ostream& operator << (std::ostream &os, const tMatBase &m)
    {
        F(j,R) {
            os << (j == 0 ? '[' : ' ');
            F(i,C) os << ' ' << m(j,i);
            os << ' ' << (j == R - 1 ? ']' : '\n');
        }
        return os << std::flush;
    }
    
    typedef tMatBlock<S, tMatBase, R, C, 1, C> Row;
    
    IL Row row(U i) { return Row(*this, i, 0); }
    
    IL const Row row(U i) const
    { return Row(const_cast<tMatBase &>(*this), i, 0); }
    
    typedef tMatBlock<S, tMatBase, R, C, R, 1> Col;
    IL Col col(U i) { return Col(*this, 0, i); }
    
    IL const Col col(U i) const
    { return Col(const_cast<tMatBase &>(*this), 0, i); }
    
    template <U NumRows, U NumCols>
    IL tMatBlock<S, tMatBase, R, C, NumRows, NumCols>
    block(U offsetRow, U offsetColumn)
    {
        return tMatBlock<S, tMatBase, R, C, NumRows, NumCols>
        (*this, offsetRow, offsetColumn);
    }
    
    template <U NumRows, U NumCols>
    IL const tMatBlock<S, tMatBase, R, C, NumRows, NumCols>
    block(U offsetRow, U offsetColumn) const
    {
        return tMatBlock<S, tMatBase, R, C, NumRows, NumCols>
        (const_cast<tMatBase &>(*this), offsetRow, offsetColumn);
    }
    
    template <U NumRows>
    IL tMatBlock<S, tMatBase, R, C, NumRows, NumRows>
    block(U offsetRow, U offsetColumn)
    {
        return tMatBlock<S, tMatBase, R, C, NumRows, NumRows>
        (*this, offsetRow, offsetColumn);
    }
    
    template <U NumRows>
    const IL tMatBlock<S, tMatBase, R, C, NumRows, NumRows>
    block(U offsetRow, U offsetColumn) const
    {
        return tMatBlock<S, tMatBase, R, C, NumRows, NumRows>
        (const_cast<tMatBase &>(*this), offsetRow, offsetColumn);
    }
    
    IL Mat<S, R, C> operator - () const
    { Mat<S, R, C> n; F(j,C) F(i,R) n(i,j) = -MI(i,j); return n; }
    
    /// Returns the transpose
    IL Mat<S, C, R> transposed() const
    { Mat<S, C, R> n; F(j,C) F(i,R) n(j,i) = MI(i,j); return n; }
    
    typedef tMatSingle<S, R, C, true> IdentityExpression;
    // Returns an identity matrix expression.
    static IL const IdentityExpression Identity()
    { return IdentityExpression::Construct(); }
    
    typedef tMatSingle<S, R, C, false> OnesExpression;
    // Returns a matrix expression with all values set to one.
    static IL const OnesExpression Ones()
    { return OnesExpression::Construct(1,1); }
    
    typedef tMatSingle<S, R, C, false> ZeroExpression;
    // Returns a matrix expression with all values set to zero.
    static IL const ZeroExpression Zero()
    { return ZeroExpression::Construct(0,0); }
};

template <class S, U R, U C, class Derived, bool IsRow>
struct tMatBase<S, R, C, Derived, IsRow, true>
{
private:
    template <class V>
    IL V _convert() const { V v; F(j, LG(R,C)) v[j] = VI(j); return v; }
public:
    DEF_MI;
    
    typedef int DFMatBase;
    
    IL Mat<S, R, C> operator - () const
    { Mat<S, R, C> n; F(j,C) F(i,R) n(i,j) = -MI(i,j); return n; }
    
    /// Returns the transpose
    IL Mat<S, C, R> transposed() const
    { Mat<S, C, R> n; F(j,C) F(i,R) n(j,i) = MI(i,j); return n; }
    
    IL S & operator [] (int i)
    { return MI(IsRow ? 0 : i, IsRow ? i : 0); }
    
    IL S operator [] (int i) const
    { return MI(IsRow ? 0 : i, IsRow ? i : 0); }
    
    template <class V>
    IL operator V() const volatile
    { return const_cast<const tMatBase*>(this)->_convert<V>(); }
    
    template <U NumRows, U NumCols>
    IL tMatBlock<S, tMatBase, R, C, NumRows, NumCols>
    block(U offsetRow, U offsetColumn)
    {
        return tMatBlock<S, tMatBase, R, C, NumRows, NumCols>
        (*this, offsetRow, offsetColumn);
    }
    
    template <U NumRows, U NumCols>
    IL const tMatBlock<S, tMatBase, R, C, NumRows, NumCols>
    block(U offsetRow, U offsetColumn) const
    {
        return tMatBlock<S, tMatBase, R, C, NumRows, NumCols>
        (const_cast<tMatBase &>(*this), offsetRow, offsetColumn);
    }
    
    template <U Dims>
    IL tMatBlock<S, tMatBase, R, C, (IsRow ? 1 : Dims), (IsRow ? Dims : 1)>
    block(U offset)
    {
        return tMatBlock<S, tMatBase, R, C, (IsRow ? 1 : Dims), (IsRow ? Dims : 1)>
        (*this, (IsRow ? 0 : offset), (IsRow ? offset : 0));
    }
    
    template <U Dims>
    const IL tMatBlock<S, tMatBase, R, C, (IsRow ? 1 : Dims), (IsRow ? Dims : 1)>
    block(U offset) const
    {
        return tMatBlock<S, tMatBase, R, C, (IsRow ? 1 : Dims), (IsRow ? Dims : 1)>
        (const_cast<tMatBase &>(*this), (IsRow ? 0 : offset), (IsRow ? offset : 0));
    }
    
    typedef tMatTilda<S, tMatBase, R, C> Tilda;
    IL Tilda operator ~() const { return Tilda(*this); }
    
    typedef tMatScalarAssigner<tMatBase, R, C, LG(R, C) - (LG(R, C) > 0), 1>
    ScalarAssigner;
    
    template <class T>
    IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value,
    ScalarAssigner>::type operator << (const T &t)
    { MI(0,0) = t; return ScalarAssigner(*this); }
    
    friend std::ostream& operator << (std::ostream &os, const tMatBase &m)
    {
        os << '[';
        F(j, LG(R, C)) os << ' ' << m[j];
        return  os << (IsRow ? " ]" : " ] T") << std::flush;
    }
    /// Returns the dot product with v.
    template <class T, class V>
    IL S dot(const tMatBase<T, R, C, V, IsRow, true> &v) const
    { S d = 0; F(i, LG(R,C)) d += VI(i) * v[i]; return d; }
    /// Returns the dot product with v.
    template <class T, class V>
    IL S dot(const tMatBase<T, C, R, V, !IsRow, true> &v) const
    { S d = 0; F(i, LG(R,C)) d += VI(i) * v[i]; return d; }
    /// Returns the squared of the norm.
    IL S normSq() const { return dot(*this); }
    /// Returns the length of v.
    IL S norm() const { return sqrt(dot(*this)); }
    /// Homogenize the vector in place.
    IL tMatBase & homogenize()
    { F(i, LG(R,C)-1) VI(i) /= VI(LG(R,C)-1); VI(LG(R,C)-1) = 1; return *this; }
    /// Returns a homogenized copy of vector.
    IL Mat<S, R, C> homogenized() const
    { Mat<S, R, C> u(*this); return u.homogenize(); }
    /// Normalizes the vector in place.
    IL tMatBase & normalize()
    { S invNorm = 1 / norm(); F(i, LG(R,C)) VI(i) *= invNorm; return *this; }
    /// Returns a normalized copy of vector.
    IL Mat<S, R, C> normalized() const
    { Mat<S, R, C> u(*this); return u.normalize(); }
    
    typedef tMatSingle<S, R, C, true> IdentityExpression;
    // Returns an identity matrix expression.
    static IL const IdentityExpression Identity()
    { return IdentityExpression::Construct(); }
    
    typedef tMatSingle<S, R, C, false> OnesExpression;
    // Returns a matrix expression with all values set to one.
    static IL const OnesExpression Ones()
    { return OnesExpression::Construct(1,1); }
    
    typedef tMatSingle<S, R, C, false> ZeroExpression;
    // Returns a matrix expression with all values set to zero.
    static IL const ZeroExpression Zero()
    { return ZeroExpression::Construct(0,0); }
};

// Identity Expressions --------------------------------------------------------

//~~~~~~~~~~~~~~~~~~~~~~~~~~
#define DEF_MAT_SINGLE(R, C) \
\
private: S _,g; IL tMatSingle(S _, S g): _(_), g(g) {}; public: \
static IL const tMatSingle Construct(S _, S g=0) { return tMatSingle(_, g); } \
IL S operator () (U row, U col) const \
{ return R==1 || C==1 ? _ : row == col ? _ : g; } \
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define DEF_MAT_SINGLE_IDENTITY(R, C) \
\
private: IL tMatSingle() {}; public: \
static IL const tMatSingle Construct() { return tMatSingle(); } \
IL S operator () (U row, U col) const \
{ return R==1 || C==1 ? 1 : row == col ? 1 : 0; } \

template <class S, U R, U C> struct tMatSingle<S, R, C, false>:
public tMatBase<S, R, C, tMatSingle<S, R, C, false>, R==1, R==1 || C==1>
{ DEF_MAT_SINGLE(R, C); };

template <class S, U R, U C> struct tMatSingle<S, R, C, true>:
public tMatBase<S, R, C, tMatSingle<S, R, C, true>, R==1, R==1 || C==1>
{ DEF_MAT_SINGLE_IDENTITY(R, C); };

template <class S, U R> struct tMatSingle<S, R, R, false>:
public tMatBase<S, R, R, tMatSingle<S, R, R, false>, R==1, R==1>
{
    DEF_MAT_SINGLE(R, R);
    /// Returns the determinant
    IL S determinant() const
    { if (_ == g) return 0; S d = _; F(i, R-1) d *= _; return d; }
    /// Returns the inverse
    IL tMatSingle inverse() const
    { return _ == g ? tMatSingle((S)(R==1)/0, (S)0/0) : tMatSingle(1/_, 0); }
};

template <class S, U R> struct tMatSingle<S, R, R, true>:
public tMatBase<S, R, R, tMatSingle<S, R, R, true>, R==1, R==1>
{
    DEF_MAT_SINGLE_IDENTITY(R, R);
    /// Returns the determinant
    IL S determinant() const { return 1; }
    /// Returns the inverse
    IL tMatSingle inverse() const { return tMatSingle(); }
};

#define DEF_VEC_SINGLE(DIMS, R, C) \
\
template <class S> struct tMatSingle<S, R, C, false>: \
public tVec##DIMS##Base<S, tMatSingle<S, R, C, false>, R, C> \
{ DEF_MAT_SINGLE(R, C); }; \
\
template <class S> struct tMatSingle<S, R, C, true>: \
public tVec##DIMS##Base<S, tMatSingle<S, R, C, true>, R, C> \
{ DEF_MAT_SINGLE_IDENTITY(R, C); };

DEF_VEC_SINGLE(2,2,1); DEF_VEC_SINGLE(2,1,2);
DEF_VEC_SINGLE(3,3,1); DEF_VEC_SINGLE(3,1,3);
DEF_VEC_SINGLE(4,4,1); DEF_VEC_SINGLE(4,1,4);

#undef DEF_VEC_SINGLE
#undef DEF_MAT_SINGLE
#undef DEF_MAT_SINGLE_IDENTITY

#define DEF_MAT_SINGLE_OP_SCALAR(OP) \
\
template <class S, U R, U C, bool II, class T> \
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, \
const tMatSingle<S, R, C, 0> >::type operator OP \
(const tMatSingle<S, R, C, II> &a, const T &t) \
{ return tMatSingle<S, R, C, 0>::Construct(a(0,0) OP t, a(1,0) OP t); } \
\
template <class S, U R, U C, bool II, class T> \
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, \
const tMatSingle<S, R, C, 0> >::type operator OP \
(const T &t, const tMatSingle<S, R, C, II> &a) \
{ return tMatSingle<S, R, C, 0>::Construct(t OP a(0,0), t OP a(1,0)); }

#define DEF_MAT_SINGLE_OP(OP) \
\
DEF_MAT_SINGLE_OP_SCALAR(OP) \
\
template <class S, U R, U C, bool II0, bool II1> \
IL const tMatSingle<S, R, C, 0> operator OP \
(const tMatSingle<S, R, C, II0> &a, const tMatSingle<S, R, C, II1> &b) \
{ return tMatSingle<S, R, C, 0>::Construct(a(0,0) OP b(0,0), a(1,0) OP b(1,0));}

DEF_MAT_SINGLE_OP(+);
DEF_MAT_SINGLE_OP(-);
DEF_MAT_SINGLE_OP_SCALAR(*);
DEF_MAT_SINGLE_OP_SCALAR(/);

#undef DEF_MAT_SINGLE_OP_SCALAR
#undef DEF_MAT_SINGLE_OP

// For better ellision with identities

template <class S, U R, class M, bool IR, bool IV>
IL const tMatBase<S, R, R, M, IR, IV> & operator *
(const tMatSingle<S, R, R, true> &a, const tMatBase<S, R, R, M, IR, IV> &b)
{ return b; }

template <class S, U R, class M, bool IR, bool IV>
IL const tMatBase<S, R, R, M, IR, IV> & operator *
(const tMatBase<S, R, R, M, IR, IV> &a, const tMatSingle<S, R, R, true> &b)
{ return a; }

template <class S, U R, class M, bool IR, bool IV>
IL tMatBase<S, R, R, M, IR, IV> & operator *=
(tMatBase<S, R, R, M, IR, IV> &a, const tMatSingle<S, R, R, true> &b)
{ return a; }

template <class S, U R, U C, bool II0, bool II1>
IL const tMatSingle<S, R, R, II0 && II1> operator *
(const tMatSingle<S, R, C, II0> &a, const tMatSingle<S, C, R, II1> &b)
{ return tMatSingle<S, R, R, II0 && II1>::Construct(a.value() * b.value()); }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define DEF_MAT_BLOCK_E_WISE_OP_EQ_SCALAR(OP, R, C) \
\
template <class T> \
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, \
tMatBlock & >::type operator OP##= (const T &t) \
{ F(i,R) F(j,C) this->operator()(i,j) OP##= t; return *this; }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define DEF_MAT_BLOCK_E_WISE_OP_EQ(OP, R, C) \
\
template <class T, class M> \
IL tMatBlock & operator OP##= \
(const tMatBase<T, R, C, M, R==1, R==1 || C==1> &m) \
{ F(i,R) F(j,C) this->operator()(i,j) OP##= m(i,j); return *this; } \
\
template <class T, class M> \
IL tMatBlock & operator OP##= \
(const tMatBase<T, C, R, M, R!=1, true> &m) \
{ F(i,R) F(j,C) this->operator()(i,j) OP##= m(j,i); return *this; } \
\
DEF_MAT_BLOCK_E_WISE_OP_EQ_SCALAR(OP, R, C)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define DEF_MAT_BLOCK_E_WISE_OP_EQ_TILDA(OP, R, C) \
\
template <class T, class M> \
IL tMatBlock & operator OP##= (const tMatTilda<T, M, R, C> &m) \
{ F(j,C) F(i,R) MI(i,j) OP##= m(i,j); return *this; }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define DEF_MAT_BLOCK_COMMON(R, C) \
\
private: Derived &_d; U _r, _c; public: \
\
template <class T, class M> \
IL tMatBlock & operator = (const tMatBase<T, LG(R,C), 1, M, false, true> &m) \
{ F(j,LG(R,C)) VI(j) = m[j]; return *this; } \
\
template <class T, class M> \
IL tMatBlock & operator = (const tMatBase<T, 1, LG(R,C), M, true, true> &m) \
{ F(j,LG(R,C)) VI(j) = m[j]; return *this; } \
\
template <class T, class M, bool IR> \
IL tMatBlock & operator = (const tMatBase<T, R, C, M, IR, false> &m) \
{ F(j,C) F(i,R) MI(i,j) = m(i,j); return *this; } \
\
template <class T> \
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, \
tMatBlock & >::type operator = (const T &t) \
{ F(j,C) F(i,R) MI(i,j) = t; return *this; } \
\
template <class T, class M> \
IL tMatBlock & operator = (const tMatBlock &m) \
{ F(j,C) F(i,R) MI(i,j) = m(i,j); return *this; } \
\
IL S & operator()(U row, U col) \
{ static S g = 0; \
return _r+row >= MR || _c+col >= MC ? g : _d(_r+row,_c+col); } \
\
IL S operator()(U row, U col) const \
{ return _r+row >= MR || _c+col >= MC ? 0 : _d(_r+row,_c+col); } \
\
DEF_MAT_BLOCK_E_WISE_OP_EQ(+, R, C) \
DEF_MAT_BLOCK_E_WISE_OP_EQ(-, R, C) \
DEF_MAT_BLOCK_E_WISE_OP_EQ_SCALAR(/, R, C) \
DEF_MAT_BLOCK_E_WISE_OP_EQ_SCALAR(*, R, C) \
DEF_MAT_BLOCK_E_WISE_OP_EQ_TILDA(/, R, C) \
DEF_MAT_BLOCK_E_WISE_OP_EQ_TILDA(*, R, C)

// Mat -------------------------------------------------------------------------

template <class S, class Derived, U MR, U MC, U R, U C> struct tMatBlock:
public tMatBase<S, R, C, tMatBlock<S, Derived, MR, MC, R, C>,
R==1, R==1 || C==1>
{
    IL tMatBlock(Derived &d, U r, U c): _d(d), _r(r), _c(c) {}
    
    DEF_MAT_BLOCK_COMMON(R, C);
};

template <class S, U NumRows, U NumCols = NumRows> struct Mat:
public tMatBase<S, NumRows, NumCols, Mat<S, NumRows, NumCols>,
NumRows == 1, NumRows == 1 || NumCols == 1>
{
    S _[NumRows * NumCols];
    
    /// Constructs a matrix with the diagonal = s
    IL Mat(S s = 0) { tMatUtil::InitMat<NumRows, NumCols>(_, s); }
    
    template <class T>
    IL Mat(const T &t) { tMatUtil::InitMat<NumRows, NumCols>(_, t); }
    
    IL Mat(const Mat &m) { F(j, NumRows * NumCols) _[j] = m._[j]; }
    
    IL operator S * () { return _; }
    IL operator const S * () const { return _; }
    
    IL S & operator()(U row, U col) { return _[col * NumRows + row]; }
    IL S operator()(U row, U col) const
    { return _[col * NumRows + row]; }
};

#define DI(IS_ROW, POS) _d(IS_ROW ? 0 : _r + POS, IS_ROW ? _c + POS : 0)

#define DEF_VEC_COMMON(CLASS, DIMS, SEL) \
IL S & operator()(U row, U col) { return _[SEL]; } \
IL S operator()(U row, U col) const { return _[SEL]; } \
IL CLASS(S s = 0) { F(j, DIMS) _[j] = s; } \
IL CLASS(const CLASS &m) { F(j, DIMS) _[j] = m._[j]; } \
template <class V> CLASS(const V &v) { tMatUtil::InitVec<DIMS>(_, v); } \
IL operator S * () { return _; } \
IL operator const S * () const { return _; }

// Vec -------------------------------------------------------------------------

template <class S, U Dims> struct Vec:
public tMatBase<S, Dims, 1, Vec<S, Dims>, false, true>
{
    S _[Dims];
    DEF_VEC_COMMON(Vec, Dims, row);
};

// Vec2 ------------------------------------------------------------------------

template <class S, class Derived, U R, U C> struct tVec2Base:
public tMatBase<S, R, C, tVec2Base<S, Derived, R, C>, R==1, true>
{
    DEF_MI;
    /// Returns the cross product with v.
    template <class T, class V>
    IL Vec<S, 3> cross(const tMatBase<T, 2, 1, V, false, true> &v) const
    {
        return Vec<S, 3>(0, 0, VI(0) * v[1] - VI(1) * v[0]);
    }
    /// Returns the cross product with v.
    template <class T, class V>
    IL Vec<S, 3> cross(const tMatBase<T, 1, 2, V, true, true> &v) const
    {
        return Vec<S, 3>(0, 0, VI(0) * v[1] - VI(1) * v[0]);
    }
    /// Returns the Up vector where y = 1 (OpenGL coordinates)
    IL static Vec<S, 2> Up() { return Vec<S, 2>(0, 1); }
    /// Returns the Right vector where x = 1 (OpenGL coordinates)
    IL static Vec<S, 2> Right() { return Vec<S, 2>(1, 0); }
};

#define DEF_VEC_BLOCK(R, C) DEF_MAT_BLOCK_COMMON(R, C); \
S &x,&y, &r,&g, &s,&t; \
IL tMatBlock(Derived &_d, U _r, U _c): _d(_d), _r(_r), _c(_c), \
x(DI(R<C,0)), y(DI(R<C,1)), \
r(DI(R<C,0)), g(DI(R<C,1)), \
s(DI(R<C,0)), t(DI(R<C,1)) {}

template <class S, class Derived, U MR, U MC>
struct tMatBlock<S, Derived, MR, MC, 2, 1>:
public tVec2Base<S, tMatBlock<S, Derived, MR, MC, 2, 1>, 2, 1>
{ DEF_VEC_BLOCK(2, 1); };

template <class S, class Derived, U MR, U MC>
struct tMatBlock<S, Derived, MR, MC, 1, 2>:
public tVec2Base<S, tMatBlock<S, Derived, MR, MC, 1, 2>, 1, 2>
{ DEF_VEC_BLOCK(1, 2); };

#undef DEF_VEC_BLOCK

#define DEF_VEC(CLASS, SEL) DEF_VEC_COMMON(CLASS, 2, SEL) \
union { S _[2]; struct { S x,y; }; struct { S r,g; }; struct { S s,t; }; }; \
IL CLASS(S x, S y): x(x), y(y) {}

template <class S> struct Mat <S, 2, 1>:
public tVec2Base<S, Mat<S, 2, 1>, 2, 1> { DEF_VEC(Mat, row); };

template <class S> struct Mat <S, 1, 2>:
public tVec2Base<S, Mat<S, 1, 2>, 1, 2> { DEF_VEC(Mat, col); };

template <class S> struct Vec <S, 2>:
public tVec2Base<S, Vec<S, 2>, 2, 1> { DEF_VEC(Vec, row); };

#undef DEF_VEC

// Vec3 ------------------------------------------------------------------------

template <class S, class Derived, U R, U C> struct tVec3Base:
public tMatBase<S, R, C, tVec3Base<S, Derived, R, C>, R==1, true>
{
    DEF_MI;
    /// Returns the cross product with v.
    template <class T, class V>
    IL Vec<S, 3> cross(const tMatBase<T, 3, 1, V, false, true> &v) const
    {
        return Vec<S, 3>(VI(1) * v[2] - VI(2) * v[1],
                         VI(2) * v[0] - VI(0) * v[2],
                         VI(0) * v[1] - VI(1) * v[0]);
    }
    /// Returns the cross product with v.
    template <class T, class V>
    IL Vec<S, 3> cross(const tMatBase<T, 1, 3, V, true, true> &v) const
    {
        return Vec<S, 3>(VI(1) * v[2] - VI(2) * v[1],
                         VI(2) * v[0] - VI(0) * v[2],
                         VI(0) * v[1] - VI(1) * v[0]);
    }
    /// Returns the Up vector where y = 1 (OpenGL coordinates)
    IL static Vec<S, 3> Up() { return Vec<S, 3>(0, 1, 0); }
    /// Returns the Right vector where x = 1 (OpenGL coordinates)
    IL static Vec<S, 3> Right() { return Vec<S, 3>(1, 0, 0); }
    /// Returns the Back vector where z = -1 (OpenGL coordinates)
    IL static Vec<S, 3> Forward() { return Vec<S, 3>(0, 0, -1); }
};

#define DEF_VEC_BLOCK(R, C) DEF_MAT_BLOCK_COMMON(R, C) \
S &x,&y,&z, &r,&g,&b, &s,&t,&p; \
IL tMatBlock(Derived &_d, U _r, U _c): _d(_d), _r(_r), _c(_c), \
x(DI(R<C,0)), y(DI(R<C,1)), z(DI(R<C,2)), \
r(DI(R<C,0)), g(DI(R<C,1)), b(DI(R<C,2)), \
s(DI(R<C,0)), t(DI(R<C,1)), p(DI(R<C,2)) {}

template <class S, class Derived, U MR, U MC>
struct tMatBlock<S, Derived, MR, MC, 3, 1>:
public tVec3Base<S, tMatBlock<S, Derived, MR, MC, 3, 1>, 3, 1>
{ DEF_VEC_BLOCK(3, 1); };

template <class S, class Derived, U MR, U MC>
struct tMatBlock<S, Derived, MR, MC, 1, 3>:
public tVec3Base<S, tMatBlock<S, Derived, MR, MC, 1, 3>, 1, 3>
{ DEF_VEC_BLOCK(1, 3); };

#undef DEF_VEC_BLOCK

#define DEF_VEC(CLASS, SEL) DEF_VEC_COMMON(CLASS, 3, SEL) union { S _[3]; \
struct { S x,y,z; }; struct { S r,g,b; }; struct { S s,t,p; }; }; \
IL CLASS(S x, S y, S z): x(x), y(y), z(z) {} \
\
template <class S0, class S1, class D, U R, U C> \
IL CLASS(const tVec2Base<S0, D, R, C> &xy, S1 z): \
x(xy[0]), y(xy[1]), z(z) {}\
\
template <class S0, class S1, class D, U R, U C> \
IL CLASS(S0 x, const tVec2Base<S1, D, R, C> &yz): \
x(x), y(yz[0]), z(yz[1]) {} \

template <class S> struct Mat <S, 3, 1>:
public tVec3Base<S, Mat<S, 3, 1>, 3, 1> { DEF_VEC(Mat, row); };

template <class S> struct Mat <S, 1, 3>:
public tVec3Base<S, Mat<S, 1, 3>, 1, 3> { DEF_VEC(Mat, col); };

template <class S> struct Vec <S, 3>:
public tVec3Base<S, Vec<S, 3>, 3, 1> { DEF_VEC(Vec, row); };

#undef DEF_VEC

// Vec4 ------------------------------------------------------------------------

template <class S, class Derived, U R, U C> struct tVec4Base:
public tMatBase<S, R, C, tVec4Base<S, Derived, R, C>, R==1, true>
{
    DEF_MI;
    /// Returns the cross product with v.
    template <class T, class V>
    IL Vec<S, 4> cross(const tMatBase<T, 4, 1, V, false, true> &v) const
    {
        return Vec<S, 4>(VI(1) * v[2] - VI(2) * v[1],
                         VI(2) * v[0] - VI(0) * v[2],
                         VI(0) * v[1] - VI(1) * v[0], 0);
    }
    /// Returns the cross product with v.
    template <class T, class V>
    IL Vec<S, 4> cross(const tMatBase<T, 1, 4, V, true, true> &v) const
    {
        return Vec<S, 4>(VI(1) * v[2] - VI(2) * v[1],
                         VI(2) * v[0] - VI(0) * v[2],
                         VI(0) * v[1] - VI(1) * v[0], 0);
    }
    /// Returns the Up vector where y = 1 (OpenGL coordinates)
    IL static Vec<S, 4> Up() { return Vec<S, 4>(0, 1, 0, 0); }
    /// Returns the Right vector where x = 1 (OpenGL coordinates)
    IL static Vec<S, 4> Right() { return Vec<S, 4>(1, 0, 0, 0); }
    /// Returns the Back vector where z = -1 (OpenGL coordinates)
    IL static Vec<S, 4> Forward() { return Vec<S, 4>(0, 0, -1, 0); }
};

#define DEF_VEC_BLOCK(R, C) DEF_MAT_BLOCK_COMMON(R, C) \
S &x,&y,&z,&w, &r,&g,&b,&a, &s,&t,&p,&q; \
IL tMatBlock(Derived &_d, U _r, U _c): _d(_d), _r(_r), _c(_c), \
x(DI(R<C,0)), y(DI(R<C,1)), z(DI(R<C,2)), w(DI(R<C,3)), \
r(DI(R<C,0)), g(DI(R<C,1)), b(DI(R<C,2)), a(DI(R<C,3)), \
s(DI(R<C,0)), t(DI(R<C,1)), p(DI(R<C,2)), q(DI(R<C,3)) {}

template <class S, class Derived, U MR, U MC>
struct tMatBlock<S, Derived, MR, MC, 4, 1>:
public tVec4Base<S, tMatBlock<S, Derived, MR, MC, 4, 1>, 4, 1>
{ DEF_VEC_BLOCK(4, 1); };

template <class S, class Derived, U MR, U MC>
struct tMatBlock<S, Derived, MR, MC, 1, 4>:
public tVec4Base<S, tMatBlock<S, Derived, MR, MC, 1, 4>, 1, 4>
{ DEF_VEC_BLOCK(1, 4); };

#undef DEF_VEC_BLOCK

#define DEF_VEC(CLASS, SEL) DEF_VEC_COMMON(CLASS, 4, SEL) union { S _[4]; \
struct { S x,y,z,w; }; struct { S r,g,b,a; }; struct { S s,t,p,q; }; }; \
IL CLASS(S x, S y, S z, S w): x(x), y(y), z(z), w(w) {} \
\
template <class S0, class S1, class S2, class D, U R, U C> \
IL CLASS(const tVec2Base<S0, D, R, C> &xy, S1 z, S2 w): \
x(xy[0]), y(xy[1]), z(z), w(w) {}\
\
template <class S0, class S1, class S2, class D, U R, U C> \
IL CLASS(S0 x, const tVec2Base<S1, D, R, C> &yz, S2 w): \
x(x), y(yz[0]), z(yz[1]), w(w) {} \
\
template <class S0, class S1, class S2, class D, U R, U C> \
IL CLASS(S0 x, S1 y, const tVec2Base<S2, D, R, C> &zw): \
x(x), y(y), z(zw[0]), w(zw[1]) {} \
\
template <class S0, class S1, class D0, \
class D1, U R0, U R1, U C0, U C1> \
IL CLASS(const tVec2Base<S0, D0, R0, C0> &xy, \
const tVec2Base<S1, D1, R1, C1> &zw): \
x(xy[0]), y(xy[1]), z(zw[0]), w(zw[1]) {} \
\
template <class S0, class S1, class D, U R, U C> \
IL CLASS(const tVec3Base<S0, D, R, C> &xyz, S1 w): \
x(xyz[0]), y(xyz[1]), z(xyz[2]), w(w) {} \
\
template <class S0, class S1, class D, U R, U C> \
IL CLASS(S0 x, const tVec3Base<S1, D, R, C> &yzw): \
x(x), y(yzw[0]), z(yzw[1]), w(yzw[2]) {}


template <class S> struct Mat <S, 4, 1>:
public tVec4Base<S, Mat<S, 4, 1>, 4, 1> { DEF_VEC(Mat, row); };

template <class S> struct Mat <S, 1, 4>:
public tVec4Base<S, Mat<S, 1, 4>, 1, 4> { DEF_VEC(Mat, col); };

template <class S> struct Vec <S, 4>:
public tVec4Base<S, Vec<S, 4>, 4, 1> { DEF_VEC(Vec, row); };

#undef DEF_VEC
#undef DI

// Vec Ops ---------------------------------------------------------------------

#define DEF_MAT_TILDA_OP(OP) \
template <class S, U R, U C, bool IR, bool IV, class V0, class V1> \
IL Mat<S, R, C> operator OP \
(const tMatBase<S, R, C, V0, IR, IV> &u, \
const tMatTilda<S, V1, R, C> &v) \
{ Mat<S, R, C> r; F(j,C) F(i,R) r(i,j) = u(i,j) OP v(i,j); return r; } \
\
template <class S, class T, U R, U C, bool IR, bool IV, class V0, class V1> \
IL tMatBase<S, R, C, V0, IR, IV> & operator OP##= \
(tMatBase<S, R, C, V0, IR, IV> &u, const tMatTilda<T, V1, R, C> &v) \
{ F(j,C) F(i,R) u(i,j) OP##= v(i,j); return u; }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define DEF_MAT_SCALAR_OP(OP) \
\
template <class S, U R, U C, class M, bool IR, bool IV, class T> \
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, Mat<S, R, C> >::type \
operator OP (const tMatBase<S, R, C, M, IR, IV> &n, const T &t) \
{ Mat<S, R, C> r; F(j,C) F(i,R) r(i,j) = n(i,j) OP t; return r; } \
\
template <class S, U R, U C, class M, bool IR, bool IV, class T> \
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, \
tMatBase<S, R, C, M, IR, IV> >::type \
operator OP##= (tMatBase<S, R, C, M, IR, IV> &n, const T &t) \
{ F(j,C) F(i,R) n(i,j) OP##= t; return n; } \
\
template <class S, U R, U C, class M, bool IR, bool IV, class T> \
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, \
Mat<S, R, C> >::type \
operator OP (const T &t, const tMatBase<S, R, C, M, IR, IV> &n) \
{ Mat<S, R, C> r; F(j,C) F(i,R) r(i,j) = t OP n(i,j); return r; }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define DEF_SIMPLE_MAT_MAT_OP(OP) \
\
template <class S, U R, U C, \
class M0, class M1, bool IR0, bool IR1, bool IV0, bool IV1> \
IL tMatBase<S, R, C, M0, IR0, IV0> operator OP##= \
(      tMatBase<S, R, C, M0, IR0, IV0> &n, \
const tMatBase<S, R, C, M1, IR1, IV1> &m) \
{ F(j,C) F(i,R) n(i,j) OP##= m(i,j); return n; } \
\
template <class S, U R, U C, \
class M0, class M1, bool IR0, bool IR1, bool IV0, bool IV1> \
IL Mat<S, R, C> operator OP \
(const tMatBase<S, R, C, M0, IR0, IV0> &n, \
const tMatBase<S, R, C, M1, IR1, IV1> &m) \
{ Mat<S, R, C> r; F(j,C) F(i,R) r(i,j) = n(i,j) OP m(i,j); return r; }

DEF_SIMPLE_MAT_MAT_OP(-);
DEF_SIMPLE_MAT_MAT_OP(+);
DEF_MAT_SCALAR_OP(-);
DEF_MAT_SCALAR_OP(+);
DEF_MAT_SCALAR_OP(*);
DEF_MAT_SCALAR_OP(/);
DEF_MAT_TILDA_OP(*);
DEF_MAT_TILDA_OP(/);

#undef DEF_MAT_TILDA_OP
#undef DEF_MAT_SCALAR_OP
#undef DEF_SIMPLE_MAT_MAT_OP

template <class S, U R, U C, U C1, class M0, class M1,
bool IR0, bool IR1, bool IV0, bool IV1>
IL Mat<S, R, C1> operator *
(const tMatBase<S, R, C, M0, IR0, IV0> &n,
 const tMatBase<S, C, C1, M1, IR1, IV1> &m)
{
    Mat<S, R, C1> r=0;
    F(j,C1) F(i,R) F(x,C) r(i,j) += n(i,x) * m(x,j);
    return r;
}

template <class S, U R, class M0, class M1,
bool IR0, bool IR1, bool IV0, bool IV1>
IL tMatBase<S, R, R, M0, IR0, IV0> & operator *=
(tMatBase<S, R, R, M0, IR0, IV0> &n,
 const tMatBase<S, R, R, M1, IR1, IV1> &m)
{
    Mat<S, R> c(n);
    F(j, R) F(i, R) F(x, R) n(i,j) += c(i,x) * m(x,j);
    return n;
}

// Quat ------------------------------------------------------------------------

template <class S> struct Quat:
public tMatBase<S, 4, 1, Quat<S>, false, true>
{
    union { S _[4]; struct { S w, x, y, z; }; };
    
    IL Quat(): w(0), x(0), y(0), z(0) {}
    
    IL Quat(S w, S x, S y, S z): w(w), x(x), y(y), z(z) {}
    
    /// Returns the identity Quaternion. [1, 0, 0, 0]
    IL static Quat Identity() { return Quat(1, 0, 0, 0); }
    
private:
    template <class V0, class V1, class V2>
    IL void _initFromBasis(const V0 &bx, const V1 &by, const V2 &bz)
    {
        S onePlusTrace = 1 + bx[0] + by[1] + bz[2];
        
        if (onePlusTrace > 1e-5) {
            S s = sqrt(onePlusTrace) * 2;
            w = 0.25 * s;
            x = (by[2] - bz[1]) / s;
            y = (bz[0] - bx[2]) / s;
            z = (bx[1] - by[0]) / s;
        } else {
            if ((bx[0] > by[1]) & (bx[0] > bz[2])) {
                S s = sqrt(1.0 + bx[0] - by[1] - bz[2]) * 2;
                w = (bz[1] - by[2]) / s;
                x = 0.25 * s;
                y = (by[0] + bx[1]) / s;
                z = (bz[0] + bx[2]) / s;
            } else if (by[1] > bz[2]) {
                S s = sqrt(1 + by[1] - bx[0] - bz[2]) * 2;
                w = (bz[0] - bx[2]) / s;
                x = (by[0] + bx[1]) / s;
                y = 0.25 * s;
                z = (bz[1] + by[2]) / s;
            } else {
                S s = sqrt(1 + bz[2] - bx[0] - by[1]) * 2;
                w = (by[0] - bx[1]) / s;
                x = (bz[0] + bx[2]) / s;
                y = (bz[1] + by[2]) / s;
                z = 0.25 * s;
            }
        }
        this->normalize();
    }
public:
    
    IL S & operator()(U row, U col) { return _[row]; }
    IL S operator()(U row, U col) const { return _[row]; }
    
    /// Constructs a quarternion from a rotation matrix.
    template <class T, class M, bool IR, bool IV>
    IL Quat(const tMatBase<T, 3, 3, M, IR, IV> &m)
    { _initFromBasis(m.col(0), m.col(1), m.col(2)); }
    
    /// Constructs a quarternion from a rotation matrix.
    template <class T, class M, bool IR, bool IV>
    IL Quat(const tMatBase<T, 4, 4, M, IR, IV> &m)
    { _initFromBasis(m.col(0), m.col(1), m.col(2)); }
    
    /// Constructs a quarternion from a rotation matrix
    template <class T0, class V0, U R0, U C0,
    class T1, class V1, U R1, U C1,
    class T2, class V2, U R2, U C2>
    IL Quat(const tVec3Base<T0, V0, R0, C0> &basisX,
            const tVec3Base<T1, V1, R1, C1> &basisY,
            const tVec3Base<T2, V2, R2, C2> &basisZ)
    { _initFromBasis(basisX, basisY, basisZ); }
    
    IL Quat & conjugate()
    { F(i, 3) _[i+1] = -_[i+1]; return *this; }
    
    IL Quat conjugated() const
    { return Quat(_[0], -_[1], -_[2], -_[3]); }
    
    IL Quat inverse() const
    { return conjugated() * (1 / this->normSq()); }
    
    IL Quat & invert()
    { return *this = inverse(); }
    
    IL Quat & log() const
    {
        S l = this->norm();
        return (Quat(0, _[1], _[2], _[3]) * (l > 1e-6 ? acos(_[0]) / l : 1));
    }
    
    IL Quat & exp() const
    {
        S t = this->norm();
        if (t < 1e-6) {
            return Quat(cos(t), _[1], _[2], _[3]);
        } else {
            S c = sin(t) / t;
            return Quat(cos(t), _[1] * c, _[2] * c, _[3] * c);
        }
    }
    
    struct AxisAngle
    {
        Vec<S, 3> axis;
        S radians;
        AxisAngle() {}
        AxisAngle(const Vec<S, 3> &axis, S radians):
        axis(axis), radians(radians) {}
    };
    
    IL AxisAngle getAxisAngle()
    { return AxisAngle(acos(w)*2, Vec<S, 3>(x,y,z).normalized()); }
    
    template <class T, class V, U R, U C>
    IL Quat & setAxisAngle(S radians, const tVec3Base<T, V, R, C> &axis)
    {
        S r = radians / 2, s = sin(r) / axis.norm();
        _[0] = cos(r); _[1] = axis.x * s; _[2] = axis.y * s; _[3] = axis.z * s;
    }
    
    IL Quat & setAxisAngle(const AxisAngle &axisAngle)
    { setAxisAngle(axisAngle.axis, axisAngle.radians); }
    
    template <class T>
    IL Quat(S radians, const Vec<T, 3> &axis) { setAxisAngle(radians, axis); }
    
    IL Quat(const AxisAngle &axisAngle) { setAxisAngle(axisAngle); }
};

/// Returns the linear interpolation between a and b
template <class S, class T>
IL Quat<S> lerp(const Quat<S> &a, const Quat<S> &b, T t)
{ return (b + t * (b - a)).normalized(); }

/// Returns the spherical interpolation between a and b
template <class S, class T>
IL Quat<S> slerp(const Quat<S> &a, const Quat<S> &b, T t,
                 bool allowFlip = true)
{
    S c1, c2, cosAngle = a.dot(b);
    if (1 - fabs(cosAngle) < 0.001) {
        c1 = 1 - t; c2 = t;
    } else {
        S angle = acos(fabs(cosAngle));
        S sinAngle = sin(angle);
        c1 = sin(angle * (1 - t)) / sinAngle;
        c2 = sin(angle * t) / sinAngle;
    }
    // Use the shortest path
    if (allowFlip && (cosAngle < 0)) c1 = -c1;
    return c1 * a + c2 * b;
}

/// Returns the spherical quadrangle interpolation between a and b
template <class S, class T>
IL Quat<S> squad(const Quat<S> &a, const Quat<S> &tanA,
                 const Quat<S> &tanB, const Quat<S> &b, T t)
{
    Quat<S> ab = slerp(a, b, t, true), tangent = slerp(tanA, tanB, t, false);
    return slerp(ab, tangent, 2 * t * (1 - t), false);
}

/// Returns the cubic interpolation between a and b
template <class S, class T>
IL Quat<S> cubicInterpolate(const Quat<S> &q0, const Quat<S> &q1,
                            const Quat<S> &q2, const Quat<S> &q3, T t)
{
    Quat<S> q0q1 = slerp(q0, q1, t + 1);
    Quat<S> q1q2 = slerp(q1, q2, t);
    Quat<S> q2q3 = slerp(q2, q3, t - 1);
    Quat<S> q0q1_q1q2 = slerp(q0q1, q1q2, 0.5 * (t + 1));
    Quat<S> q1q2_q2q3 = slerp(q1q2, q2q3, 0.5 * t);
    return slerp(q0q1_q1q2, q1q2_q2q3, t);
}

/// Returns the log difference between a and b
template <class S>
IL Quat<S> logDifference(const Quat<S> &a, const Quat<S> &b)
{ return (a.inverse() * b).normalized().log(); }

/// Returns the tangent for spherical quadrangle interpolation
template <class S>
IL Quat<S> squadTangent(const Quat<S> &before,
                        const Quat<S> &center,
                        const Quat<S> &after)
{
    Quat<S> l1 = logDifference(center, before);
    Quat<S> l2 = logDifference(center, after);
    return center * (-0.25 * (l1 + l2)).exp();
}

template <class S> IL Quat<S> operator * (const Quat<S> &p, const Quat<S> &q)
{
    return Quat<S>(p.w*q.w - p.x*q.x - p.y*q.y - p.z*q.z,
                   p.w*q.x + p.x*q.w + p.y*q.z - p.z*q.y,
                   p.w*q.y - p.x*q.z + p.y*q.w + p.z*q.x,
                   p.w*q.z + p.x*q.y - p.y*q.x + p.z*q.w);
}

template <class S> IL Quat<S> & operator *= (Quat<S> &p, const Quat<S> &q)
{ return p = p * q; }

namespace tMatUtil
{
    template <class S, U Dims> struct MatInv;
    template <class S, U Dims> struct MatDet;
    
#define SM(I0,I1,I2) (src[0x##I0] * src[0x##I1] * src[0x##I2])
#define SI(I0,I1,I2,I3,I4,I5,I6,I7,I8,I9,IA,IB,IC,ID,IE,IF,IG,IH) \
(SM(I0,I1,I2) - SM(I3,I4,I5) + SM(I6,I7,I8) - \
SM(I9,IA,IB) + SM(IC,ID,IE) - SM(IF,IG,IH))
    
    template <class S> struct MatInv <S, 4>
    {
        IL void operator()(const S *src, S *dst)
        {
            S dst0  = SI(5,a,f,5,b,e,9,7,e,9,6,f,d,6,b,d,7,a);
            S dst4  = SI(4,b,e,4,a,f,8,6,f,8,7,e,c,7,a,c,6,b);
            S dst8  = SI(4,9,f,4,b,d,8,7,d,8,5,f,c,5,b,c,7,9);
            S dst12 = SI(4,a,d,4,9,e,8,5,e,8,6,d,c,6,9,c,5,a);
            S det = src[0]*dst0 + src[1]*dst4 + src[2]*dst8 + src[3]*dst12;
            if (det == 0) { F(i, 16) dst[i] /= 0; return; }
            S invDet = 1 / det;
            dst[0] = dst0, dst[4] = dst4, dst[8] = dst8, dst[12] = dst12;
            dst[1 ] = SI(1,b,e,1,a,f,9,2,f,9,3,e,d,3,a,d,2,b);
            dst[5 ] = SI(0,a,f,0,b,e,8,3,e,8,2,f,c,2,b,c,3,a);
            dst[9 ] = SI(0,b,d,0,9,f,8,1,f,8,3,d,c,3,9,c,1,b);
            dst[13] = SI(0,9,e,0,a,d,8,2,d,8,1,e,c,1,a,c,2,9);
            dst[2 ] = SI(1,6,f,1,7,e,5,3,e,5,2,f,d,2,7,d,3,6);
            dst[6 ] = SI(0,7,e,0,6,f,4,2,f,4,3,e,c,3,6,c,2,7);
            dst[10] = SI(0,5,f,0,7,d,4,3,d,4,1,f,c,1,7,c,3,5);
            dst[14] = SI(0,6,d,0,5,e,4,1,e,4,2,d,c,2,5,c,1,6);
            dst[3 ] = SI(1,7,a,1,6,b,5,2,b,5,3,a,9,3,6,9,2,7);
            dst[7 ] = SI(0,6,b,0,7,a,4,3,a,4,2,b,8,2,7,8,3,6);
            dst[11] = SI(0,7,9,0,5,b,4,1,b,4,3,9,8,3,5,8,1,7);
            dst[15] = SI(0,5,a,0,6,9,4,2,9,4,1,a,8,1,6,8,2,5);
            F(i, 16) dst[i] *= invDet;
        }
    };
    
    template <class S> struct MatDet <S, 4>
    {
        IL S operator()(const S *src)
        {
            return (src[0] * SI(5,a,f,5,b,e,9,7,e,9,6,f,d,6,b,d,7,a) +
                    src[1] * SI(4,b,e,4,a,f,8,6,f,8,7,e,c,7,a,c,6,b) +
                    src[2] * SI(4,9,f,4,b,d,8,7,d,8,5,f,c,5,b,c,7,9) +
                    src[3] * SI(4,a,d,4,9,e,8,5,e,8,6,d,c,6,9,c,5,a) );
        }
    };
    
#undef SM
#undef SI
    
    // Mat 4x4 SSE inverse
    
#if defined(__SSE2__) || defined(_M_IX86_FP)
    // Eigen's 4x4 float inverse
#define SM _mm_mul_ps
#define SI(V, OP) _mm_shuffle_ps(V, V, OP)
#define DEF_SSE_INV \
__m128 L1,L2,L3,L4,A,B,C,D,iA,iB,iC,iD,DC,AB,dA,dB,dC,dD,dt,d,d1,d2,rd;L1=_mm_l\
oadu_ps(src+0);L2=_mm_loadu_ps(src+4);L3=_mm_loadu_ps(src+8);L4=_mm_loadu_ps(sr\
c+12);A=_mm_unpacklo_ps(L1,L2);B=_mm_unpacklo_ps(L3,L4);C=_mm_unpackhi_ps(L1,L2\
);D=_mm_unpackhi_ps(L3,L4);AB=SM(SI(A,0x0F),B);AB=_mm_sub_ps(AB,SM(SI(A,0xA5),S\
I(B,0x4E)));DC=SM(SI(D,0x0F),C);DC=_mm_sub_ps(DC,SM(SI(D,0xA5),SI(C,0x4E)));dA=\
SM(SI(A,0x5F),A);dA=_mm_sub_ss(dA,_mm_movehl_ps(dA,dA));dB=SM(SI(B,0x5F),B);dB=\
_mm_sub_ss(dB,_mm_movehl_ps(dB,dB));dC=SM(SI(C,0x5F),C);dC=_mm_sub_ss(dC,_mm_mo\
vehl_ps(dC,dC));dD=SM(SI(D,0x5F),D);dD=_mm_sub_ss(dD,_mm_movehl_ps(dD,dD));d=SM\
(SI(DC,0xD8),AB);iD=SM(SI(C,0xA0),_mm_movelh_ps(AB,AB));iD=_mm_add_ps(iD,SM(SI(\
C,0xF5),_mm_movehl_ps(AB,AB)));iA=SM(SI(B,0xA0),_mm_movelh_ps(DC,DC));iA=_mm_ad\
d_ps(iA,SM(SI(B,0xF5),_mm_movehl_ps(DC,DC)));d=_mm_add_ps(d,_mm_movehl_ps(d,d))\
;d=_mm_add_ss(d,SI(d,1));d1=_mm_mul_ss(dA,dD);d2=_mm_mul_ss(dB,dC);iD=_mm_sub_p\
s(SM(D,SI(dA,0)),iD);iA=_mm_sub_ps(SM(A,SI(dD,0)),iA);dt=_mm_sub_ss(_mm_add_ss(\
d1,d2),d);rd=_mm_div_ss(_mm_set_ss(1.0f),dt);iB=SM(D,SI(AB,0x33));iB=_mm_sub_ps\
(iB,SM(SI(D,0xB1),SI(AB,0x66)));iC=SM(A,SI(DC,0x33));iC=_mm_sub_ps(iC,SM(SI(A,0\
xB1),SI(DC,0x66)));rd=SI(rd,0);rd=_mm_xor_ps(rd,_mm_castsi128_ps(_mm_set_epi32(\
0,1<<31,1<<31,0)));iB=_mm_sub_ps(SM(C,SI(dB,0)),iB);iC=_mm_sub_ps(SM(B,SI(dC,0)\
),iC);iA=SM(rd,iA);iB=SM(rd,iB);iC=SM(rd,iC);iD=SM(rd,iD);_mm_storeu_ps(dst+0,S\
I(iB,0x77));_mm_storeu_ps(dst+4,SI(iB,0x22));_mm_storeu_ps(dst+8,SI(iD,0x77));_\
mm_storeu_ps(dst+12,SI(iD,0x22));
    
    template <> struct MatInv <float, 4>
    { IL void operator()(const float *src, float *dst) { DEF_SSE_INV; } };
    
#undef DEF_SSE_INV
#undef SI
#undef SM
    
    // Eigen's 4x4 double inverse
#define SM _mm_mul_pd
#define SS _mm_sub_pd
#define SI _mm_shuffle_pd
#define ST(O,A,K,D) _mm_storeu_pd(i+O,SM(SI(i##A##2,i##A##1,K),d##D))
#define SX(A,B) _mm_xor_pd(rd,_mm_castsi128_pd(_mm_set_epi32(1<<A,0,1<<B,0)))
#define DEF_SSE_INV \
__m128d A1,A2,B1,B2,C1,C2,D1,D2,iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,DC1,DC2,AB1,AB2\
,dA,dB,dC,dD,det,d1,d2,rd;double*i=dst;A1=_mm_loadu_pd(src+0);B1=_mm_loadu_pd(s\
rc+2);A2=_mm_loadu_pd(src+4);B2=_mm_loadu_pd(src+6);C1=_mm_loadu_pd(src+8);D1=_\
mm_load_pd(src+10);C2=_mm_loadu_pd(src+12);D2=_mm_loadu_pd(src+14);dA=SI(A2,A2,\
1);dA=SM(A1,dA);dA=_mm_sub_sd(dA,SI(dA,dA,3));dB=SI(B2,B2,1);dB=SM(B1,dB);dB=_m\
m_sub_sd(dB,SI(dB,dB,3));AB1=SM(B1,SI(A2,A2,3));AB2=SM(B2,SI(A1,A1,0));AB1=SS(A\
B1,SM(B2,SI(A1,A1,3)));AB2=SS(AB2,SM(B1,SI(A2,A2,0)));dC=SI(C2,C2,1);dC=SM(C1,d\
C);dC=_mm_sub_sd(dC,SI(dC,dC,3));dD=SI(D2,D2,1);dD=SM(D1,dD);dD=_mm_sub_sd(dD,S\
I(dD,dD,3));DC1=SM(C1,SI(D2,D2,3));DC2=SM(C2,SI(D1,D1,0));DC1=SS(DC1,SM(C2,SI(D\
1,D1,3)));DC2=SS(DC2,SM(C1,SI(D2,D2,0)));d1=SM(AB1,SI(DC1,DC2,0));d2=SM(AB2,SI(\
DC1,DC2,3));rd=_mm_add_pd(d1,d2);rd=_mm_add_sd(rd,SI(rd,rd,3));iD1=SM(AB1,SI(C1\
,C1,0));iD2=SM(AB1,SI(C2,C2,0));iD1=_mm_add_pd(iD1,SM(AB2,SI(C1,C1,3)));iD2=_mm\
_add_pd(iD2,SM(AB2,SI(C2,C2,3)));iA1=SM(DC1,SI(B1,B1,0));iA2=SM(DC1,SI(B2,B2,0)\
);iA1=_mm_add_pd(iA1,SM(DC2,SI(B1,B1,3)));iA2=_mm_add_pd(iA2,SM(DC2,SI(B2,B2,3)\
));dA=SI(dA,dA,0);iD1=SS(SM(D1,dA),iD1);iD2=SS(SM(D2,dA),iD2);dD=SI(dD,dD,0);iA\
1=SS(SM(A1,dD),iA1);iA2=SS(SM(A2,dD),iA2);d1=_mm_mul_sd(dA,dD);d2=_mm_mul_sd(dB\
,dC);iB1=SM(D1,SI(AB2,AB1,1));iB2=SM(D2,SI(AB2,AB1,1));iB1=SS(iB1,SM(SI(D1,D1,1\
),SI(AB2,AB1,2)));iB2=SS(iB2,SM(SI(D2,D2,1),SI(AB2,AB1,2)));det=_mm_add_sd(d1,d\
2);det=_mm_sub_sd(det,rd);iC1=SM(A1,SI(DC2,DC1,1));iC2=SM(A2,SI(DC2,DC1,1));iC1\
=SS(iC1,SM(SI(A1,A1,1),SI(DC2,DC1,2)));iC2=SS(iC2,SM(SI(A2,A2,1),SI(DC2,DC1,2))\
);rd=_mm_div_sd(_mm_set_sd(1.0),det);rd=SI(rd,rd,0);dB=SI(dB,dB,0);iB1=SS(SM(C1\
,dB),iB1);iB2=SS(SM(C2,dB),iB2);d1=SX(31,0);d2=SX(0,31);dC=SI(dC,dC,0);iC1=SS(S\
M(B1,dC),iC1);iC2=SS(SM(B2,dC),iC2);ST(0,A,3,1);ST(4,A,0,2);ST(2,B,3,1);ST(6,B,\
0,2);ST(8,C,3,1);ST(12,C,0,2);ST(10,D,3,1);ST(14,D,0,2);
    
    template <> struct MatInv <double, 4>
    { IL void operator()(const double *src, double *dst) { DEF_SSE_INV; }; };
    
#undef DEF_SSE_INV
#undef SM
#undef SS
#undef ST
#undef SI
#undef SX
    
#endif
    
    template <class S> struct MatInv <S, 3>
    {
        IL void operator()(const S *src, S *dst)
        {
            dst[0] = + src[4] * src[8] - src[5] * src[7];
            dst[1] = - src[1] * src[8] + src[2] * src[7];
            dst[2] = + src[1] * src[5] - src[2] * src[4];
            dst[3] = - src[3] * src[8] + src[5] * src[6];
            dst[4] = + src[0] * src[8] - src[2] * src[6];
            dst[5] = - src[0] * src[5] + src[2] * src[3];
            dst[6] = + src[3] * src[7] - src[4] * src[6];
            dst[7] = - src[0] * src[7] + src[1] * src[6];
            dst[8] = + src[0] * src[4] - src[1] * src[3];
            S det = src[0] * dst[0] + src[1] * dst[3] + src[2] * dst[6];
            det = 1 / det;
            F(i, 9) dst[i] *= det;
        }
    };
    
    template <class S> struct MatDet <S, 3>
    {
        IL S operator()(const S *src)
        {
            S d0 = + src[4] * src[8] - src[5] * src[7];
            S d3 = - src[3] * src[8] + src[5] * src[6];
            S d6 = + src[3] * src[7] - src[4] * src[6];
            return src[0] * d0 + src[1] * d3 + src[2] * d6;
        }
    };
    
    template <class S> struct MatInv<S, 2>
    {
        IL void operator()(const S *src, S *dst)
        {
            dst[0] = + src[3]; dst[1] = - src[1];
            dst[2] = - src[2]; dst[3] = + src[0];
            S det = 1 / (src[0] * dst[0] + src[1] * dst[2]);
            dst[0] *= det; dst[1] *= det;
            dst[2] *= det; dst[3] *= det;
        }
    };
    
    template <class S> struct MatDet<S, 2>
    {
        IL S operator()(const S *src)
        {
            return 1 / (src[0] * src[3] - src[2] * src[1]);
        }
    };
    
    template <class S> struct MatInv<S, 1>
    {
        IL void operator()(const S *src, S *dst) { dst[0] = 1 / src[0]; }
    };
    
    template <class S> struct MatDet<S, 1>
    {
        IL S operator()(const S *src) { return src[0]; }
    };
}

// Square Mat ------------------------------------------------------------------

//~~~~~~~~~~~~~~~~~~~~~~~~~
#define DEF_SQ_MAT_COMMON(R) \
\
template <class T> \
IL Mat(const T &t) { tMatUtil::InitMat<R, R>(_, t); } \
\
IL Mat(const Mat &m) { F(j, R*R) _[j] = m._[j]; } \
\
IL S & operator()(U row, U col) { return _[col * R + row]; } \
IL S operator()(U row, U col) const { return _[col * R + row]; } \
\
IL operator S * () { return _; } \
IL operator const S * () const { return _; } \
\
S _[R*R];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define DEF_SQ_MAT_BLOCK_COMMON(R) \
\
DEF_MAT_BLOCK_COMMON(R, R) \
\
template <class T, class M, bool IR, bool IV> \
IL tMatBlock & operator *= (const tMatBase<T, R, R, M, IR, IV> &m) \
{ Mat<S, R> c(*this); F(j, R) F(i, R) F(x, R) MI(i,j) += c(i,x) * m(x,j); \
return *this; } \
\
template <class T> \
IL tMatBlock & operator *= (const tMatSingle<T, R, R, true> &m) \
{ return *this; }

template <class S, class Derived, U R> struct tSqMatBase:
public tMatBase<S, R, R, tSqMatBase<S, Derived, R>, R==1, R==1>
{
    DEF_MI;
    
    /// Transposes the matrix or block in place
    IL tSqMatBase & transpose()
    {
        Mat<S, R, R> n(*this);
        F(i,R) F(j,R) MI(i,j) = n(j,i);
        return *this;
    }
};

template <class S, class Derived, U R> struct tInvertibleSqMatBlock:
public tSqMatBase<S, tInvertibleSqMatBlock<S, Derived, R>, R>
{
    DEF_MI;
    
    /// Invert in place
    IL tInvertibleSqMatBlock & invert()
    {
        S src[R*R], dst[R*R];
        F(j, R) F(i, R) src[j*R+i] = MI(i,j);
        tMatUtil::MatInv<S, R>()(src, dst);
        F(j, R) F(i, R) MI(i,j) = dst[j*R+i];
        return *this;
    }
    
    /// Returns the inverse
    IL Mat<S, R> inverse() const
    {
        Mat<S, R> m;
        S src[R*R]; S *dst = m;
        F(j, R) F(i, R) src[j*R+i] = MI(i,j);
        tMatUtil::MatInv<S, R>()(src, dst);
        return m;
    }
    
    /// Returns the determinant
    IL S determinant() const
    {
        S src[R*R];
        F(j, R) F(i, R) src[j*R+i] = MI(i,j);
        return tMatUtil::MatDet<S, R>()(src);
    }
};

template <class S, class Derived, U R> struct tInvertibleSqMat:
public tSqMatBase<S, tInvertibleSqMat<S, Derived, R>, R>
{
    DEF_MI;
    
    /// Invert in place
    IL tInvertibleSqMat & invert()
    {
        Mat<S, R> c(*this);
        S *dst = & MI(0,0), *src = c;
        tMatUtil::MatInv<S, R>()(src, dst);
        return *this;
    }
    
    /// Returns the inverse
    IL Mat<S, R> inverse() const
    {
        Mat<S, R> m;
        S src[R*R]; S *dst = m;
        F(j, R) F(i, R) src[j*R+i] = MI(i,j);
        tMatUtil::MatInv<S, R>()(src, dst);
        return m;
    }
    
    /// Returns the determinant
    IL S determinant() const
    {
        const S *src = & MI(0,0);
        return tMatUtil::MatDet<S, R>()(src);
    }
};

template <class S, U NumRows> struct Mat <S, NumRows, NumRows>:
public tSqMatBase<S, Mat<S, NumRows, NumRows>, NumRows>
{
    /// Constructs a matrix with the diagonal = s
    IL Mat(S s = 0)
    {
        enum { R = NumRows, RR = R*R };
        F(j, RR) _[j] = 0;
        F(j, R) _[j*R+j] = s;
    }
    
    DEF_SQ_MAT_COMMON(NumRows);
};

template <class S, class Derived, U MR, U MC, U R>
struct tMatBlock<S, Derived, MR, MC, R, R>:
public tSqMatBase<S, tMatBlock<S, Derived, MR, MC, R, R>, R>
{
    DEF_SQ_MAT_BLOCK_COMMON(R);
    IL tMatBlock(Derived &d, U r, U c): _d(d), _r(r), _c(c) {}
};

// Mat4x4 ----------------------------------------------------------------------

template <class S, class Derived, U MR, U MC>
struct tMatBlock<S, Derived, MR, MC, 4, 4>:
public tInvertibleSqMatBlock<S, tMatBlock<S, Derived, MR, MC, 4, 4>, 4>
{
    DEF_SQ_MAT_BLOCK_COMMON(4);
    IL tMatBlock(Derived &_d, U _r, U _c): _d(_d), _r(_r), _c(_c) {}
};

template <class S> struct Mat<S, 4, 4>:
public tInvertibleSqMat<S, Mat<S, 4, 4>, 4>
{
    IL Mat(S m00, S m01, S m02, S m03,
           S m10, S m11, S m12, S m13,
           S m20, S m21, S m22, S m23,
           S m30, S m31, S m32, S m33,
           bool rowWise = true)
    {
        if (rowWise) {
            _[0] = m00; _[4] = m01; _[8]  = m02; _[12] = m03;
            _[1] = m10; _[5] = m11; _[9]  = m12; _[13] = m13;
            _[2] = m20; _[6] = m21; _[10] = m22; _[14] = m23;
            _[3] = m30; _[7] = m31; _[11] = m32; _[15] = m33;
        } else {
            _[0] = m00; _[4] = m10; _[8]  = m20; _[12] = m30;
            _[1] = m01; _[5] = m11; _[9]  = m21; _[13] = m31;
            _[2] = m02; _[6] = m12; _[10] = m22; _[14] = m32;
            _[3] = m03; _[7] = m13; _[11] = m23; _[15] = m33;
        }
    }
    
    template <class V0, class V1, class V2, class V3>
    IL Mat(const tVec4Base<S, V0, 4, 1> &v0,
           const tVec4Base<S, V1, 4, 1> &v1,
           const tVec4Base<S, V2, 4, 1> &v2,
           const tVec4Base<S, V3, 4, 1> &v3) // by columns
    {
        _[0] = v0(0,0); _[4] = v1(0,0); _[8]  = v2(0,0); _[12] = v3(0,0);
        _[1] = v0(1,0); _[5] = v1(1,0); _[9]  = v2(1,0); _[13] = v3(1,0);
        _[2] = v0(2,0); _[6] = v1(2,0); _[10] = v2(2,0); _[14] = v3(2,0);
        _[3] = v0(3,0); _[7] = v1(3,0); _[11] = v2(3,0); _[15] = v3(3,0);
    }
    
    template <class V0, class V1, class V2, class V3>
    IL Mat(const tVec4Base<S, V0, 1, 4> &v0,
           const tVec4Base<S, V1, 1, 4> &v1,
           const tVec4Base<S, V2, 1, 4> &v2,
           const tVec4Base<S, V3, 1, 4> &v3) // by rows
    {
        _[0] = v0(0,0); _[4] = v0(0,1); _[8]  = v0(0,2); _[12] = v0(0,3);
        _[1] = v1(0,0); _[5] = v1(0,1); _[9]  = v1(0,2); _[13] = v1(0,3);
        _[2] = v2(0,0); _[6] = v2(0,1); _[10] = v2(0,2); _[14] = v2(0,3);
        _[3] = v3(0,0); _[7] = v3(0,1); _[11] = v3(0,2); _[15] = v3(0,3);
    }
    
    /// Constructs a matrix with the diagonal = s
    IL Mat(S s = 0) { F(j, 4*4) _[j] = 0; F(j, 4) _[j*4+j] = s; }
    
    DEF_SQ_MAT_COMMON(4);
};

// Mat3x3 ----------------------------------------------------------------------

template <class S, class Derived, U MR, U MC>
struct tMatBlock<S, Derived, MR, MC, 3, 3>:
public tInvertibleSqMatBlock<S, tMatBlock<S, Derived, MR, MC, 3, 3>, 3>
{
    DEF_SQ_MAT_BLOCK_COMMON(3);
    IL tMatBlock(Derived &_d, U _r, U _c): _d(_d), _r(_r), _c(_c) {}
};

template <class S> struct Mat<S, 3, 3>:
public tInvertibleSqMat<S, Mat<S, 3, 3>, 3>
{
    IL Mat(S m00, S m01, S m02,
           S m10, S m11, S m12,
           S m20, S m21, S m22,
           bool rowWise = true)
    {
        if (rowWise) {
            _[0] = m00; _[3] = m01; _[6] = m02;
            _[1] = m10; _[4] = m11; _[7] = m12;
            _[2] = m20; _[5] = m21; _[8] = m22;
        } else {
            _[0] = m00; _[3] = m10; _[6] = m20;
            _[1] = m01; _[4] = m11; _[7] = m21;
            _[2] = m02; _[5] = m12; _[8] = m22;
        }
    }
    
    /// Constructs a matrix with the diagonal = s
    IL Mat(S s = 0) { F(j, 3*3) _[j] = 0; F(j, 3) _[j*3+j] = s; }
    
    DEF_SQ_MAT_COMMON(3);
};

// Mat2x2 ----------------------------------------------------------------------

template <class S, class Derived, U MR, U MC>
struct tMatBlock<S, Derived, MR, MC, 2, 2>:
public tInvertibleSqMatBlock<S, tMatBlock<S, Derived, MR, MC, 2, 2>, 2>
{
    DEF_SQ_MAT_BLOCK_COMMON(2);
    IL tMatBlock(Derived &_d, U _r, U _c): _d(_d), _r(_r), _c(_c) {}
};

template <class S> struct Mat<S, 2, 2>:
public tInvertibleSqMat<S, Mat<S, 2, 2>, 2>
{
    IL Mat(S m00, S m01,
           S m10, S m11,
           bool rowWise = true)
    {
        if (rowWise) {
            _[0] = m00; _[2] = m01;
            _[1] = m10; _[3] = m11;
        } else {
            _[0] = m00; _[2] = m10;
            _[1] = m01; _[3] = m11;
        }
    }
    
    /// Constructs a matrix with the diagonal = s
    IL Mat(S s = 0) { F(j, 2*2) _[j] = 0; F(j, 2) _[j*2+j] = s; }
    
    DEF_SQ_MAT_COMMON(2);
};

// The single vector or matrix -------------------------------------------------

#define DEF_VEC \
IL S & operator()(U row, U col) { return _; } \
IL S operator()(U row, U col) const { return _; } \
IL operator S & () { return _; } \
IL operator const S & () const { return _; } \
IL operator S * () { return &_; } \
IL operator const S * () const { return &_; }

template <class S> struct Vec<S, 1>:
public tInvertibleSqMat<S, Vec<S, 1>, 1>
{
    S _;
    IL Vec(S s = 0): _(s) {};
    template <class V> Vec(const V &v) { tMatUtil::InitSingle(_, v); };
    DEF_VEC;
};

template <class S, class Derived, U MR, U MC>
struct tMatBlock<S, Derived, MR, MC, 1, 1>:
public tMatBase<S, 1, 1, tMatBlock<S, Derived, MR, MC, 1, 1>, true, true>
{
    DEF_SQ_MAT_BLOCK_COMMON(1);
    IL tMatBlock(Derived &_d, U _r, U _c): _d(_d), _r(_r), _c(_c) {}
};

template <class S> struct Mat<S, 1, 1>:
public tInvertibleSqMat<S, Vec<S, 1>, 1>
{
    S _;
    IL Mat(S s = 0): _(s) {};
    template <class V> Mat(const V &v) { tMatUtil::InitSingle(_, v); }
    DEF_VEC;
};

// Utility ----------------------------------------------------------------------


/// Rotates a 2x2 transform matrix counter-clockwise.
template <class S, class M, class A>
IL Mat<S, 2> rotate(const tMatBase<S, 2, 2, M, false, false> &m, A radians)
{
    const S c = cos(radians), s = sin(radians);
    return Mat<S, 2>(c, -s, s, c) * m;
}

/// Rotates a 3x3 transform matrix counter-clockwise about the Y axis.
template <class S, class M, class A>
IL Mat<S, 3> rotateY(const tMatBase<S, 3, 3, M, false, false> &m, A radians)
{
    const S c = cos(radians), s = sin(radians);
    return Mat<S, 3>(c, 0, s,
                     0, 1, 0,
                     -s, 0, c) * m;
}

/// Rotates a 3x3 transform matrix counter-clockwise about the Z axis.
template <class S, class M, class A>
IL Mat<S, 3> rotateZ(const tMatBase<S, 3, 3, M, false, false> &m, A radians)
{
    const S c = cos(radians), s = sin(radians);
    return Mat<S, 3>(c, -s, 0,
                     s, c, 0,
                     0, 0, 1) * m;
}

/// Rotates a 3x3 transform matrix counter-clockwise about the specified axis.
template <class S, class M, class A, class T, class V, U R, U C>
IL Mat<S, 3> rotate(const tMatBase<S, 3, 3, M, false, false> &m, A radians,
                    const tVec3Base<T, V, R, C> &axis)
{
    const Vec<S, 3> a = axis.normalized();
    const S c = cos(radians), s = sin(radians), d = 1 - c;
    const S &x = a[0], &y = a[1], &z = a[2];
    return Mat<S, 3>(x*x*d + c,   y*x*d - z*s, z*x*d + y*s,
                     x*y*d + z*s, y*y*d + c,   z*y*d - x*s,
                     x*z*d - y*s, y*z*d + x*s, z*z*d + c) * m;
}

/// Scales a 3x3 transform matrix.
template <class S, class M, class T, class V, U R, U C>
IL Mat<S, 3> scale(const tMatBase<S, 3, 3, M, false, false> &m,
                   const tVec3Base<T, V, R, C> &s)
{
    return Mat<S, 3>(s[0], 0, 0, 0,
                     0, s[1], 0, 0,
                     0, 0, s[2], 0,
                     0,  0, 0, 1) * m;
}

/// Rotates a 4x4 transform matrix counter-clockwise about the X axis.
template <class S, class M, class A>
IL Mat<S, 4> rotateX(const tMatBase<S, 4, 4, M, false, false> &m, A radians)
{
    const S c = cos(radians), s = sin(radians);
    return Mat<S, 4>(1, 0, 0, 0,
                     0, c, -s, 0,
                     0, s, c, 0,
                     0, 0, 0, 1) * m;
}

/// Rotates a 4x4 transform matrix counter-clockwise about the Y axis.
template <class S, class M, class A>
IL Mat<S, 4> rotateY(const tMatBase<S, 4, 4, M, false, false> &m, A radians)
{
    const S c = cos(radians), s = sin(radians);
    return Mat<S, 4>(c, 0, s, 0,
                     0, 1, 0, 0,
                     -s, 0, c, 0,
                     0, 0, 0, 1) * m;
}

/// Rotates a 4x4 transform matrix counter-clockwise about the Z axis.
template <class S, class M, class A>
IL Mat<S, 4> rotateZ(const tMatBase<S, 4, 4, M, false, false> &m, A radians)
{
    const S c = cos(radians), s = sin(radians);
    return Mat<S, 4>(c, -s, 0, 0,
                     s, c, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1) * m;
}

/// Rotates a 4x4 transform matrix counter-clockwise about the specified axis.
template <class S, class M, class A, class T, class V, U R, U C>
IL Mat<S, 4> rotate(const tMatBase<S, 4, 4, M, false, false> &m,
                    A radians,
                    const tVec3Base<T, V, R, C> &axis)
{
    const Vec<S, 3> a = axis.normalized();
    const S c = cos(radians), s = sin(radians), d = 1 - c;
    const S &x = a[0], &y = a[1], &z = a[2];
    return Mat<S, 4>(x*x*d + c,   y*x*d - z*s, z*x*d + y*s, 0,
                     x*y*d + z*s, y*y*d + c,   z*y*d - x*s, 0,
                     x*z*d - y*s, y*z*d + x*s, z*z*d + c,   0,
                     0, 0, 0, 1) * m;
}

/// Translates a 4x4 transform matrix.
template <class S, class M, class T, class V, U R, U C>
IL Mat<S, 4> translate(const tMatBase<S, 4, 4, M, false, false> &m,
                       const tVec3Base<T, V, R, C> &v)
{
    return Mat<S, 4>(1, 0, 0, v[0],
                     0, 1, 0, v[1],
                     0, 0, 1, v[2],
                     0, 0, 0, 1) * m;
}


/// Scales a 4x4 transform matrix.
template <class S, class M, class T, class V, U R, U C>
IL Mat<S, 4> scale(const tMatBase<S, 4, 4, M, false, false> &m,
                   const tVec3Base<T, V, R, C> &s)
{
    return Mat<S, 4>(s[0], 0, 0, 0,
                     0, s[1], 0, 0,
                     0, 0, s[2], 0,
                     0,  0, 0, 1) * m;
}

/// Constructs a look-at view matrix.
template <class S,
class V0, U R0, U C0,
class V1, U R1, U C1,
class V2, U R2, U C2>
IL Mat<S, 4> lookAt(const tVec3Base<S, V0, R0, C0> &eye,
                    const tVec3Base<S, V1, R1, C1> &center,
                    const tVec3Base<S, V2, R2, C2> &up)
{
    const Vec<S, 3> z = (eye - center).normalized();
    const Vec<S, 3> y = up;
    const Vec<S, 3> x = y.cross(z);
    return Mat<S, 4>(x[0], x[1], x[2], -x.dot(eye),
                     y[0], y[1], y[2], -y.dot(eye),
                     z[0], z[1], z[2], -z.dot(eye),
                     0, 0, 0, 1);
}

/// Constructs a orthographic view matrix.
template <class S>
IL Mat<S, 4> ortho(S width, S height, S zNear, S zFar,
                   bool directX = false)
{
    S n = zNear, f = zFar, nf = n - f;
    bool d = directX;
    return Mat<S, 4>(2 / width, 0, 0, -1,
                     0, 2 / height, 0, -1,
                     0, 0, (d ? 1 : 2) / nf, (n + (d ? 0 : f)) / nf,
                     0, 0, 0, 1);
}

/// Constructs a orthographic view matrix.
template <class S>
IL Mat<S, 4> ortho(S left, S right, S bottom, S top, S zNear, S zFar,
                   bool directX = false)
{
    S l = left, r = right, b = bottom, t = top, n = zNear, f = zFar;
    bool d = directX;
    return Mat<S, 4>(2 / (r - l), 0, 0, (l + r) / (l - r),
                     0, 2 / (t - b), 0, (t + b) / (b - t),
                     0, 0, (d ? 1 : 2) / (n - f), (n + (d ? 0 : f)) / (n - f),
                     0, 0, 0, 1);
}

/// Constructs a perspective view matrix.
template <class S>
IL Mat<S, 4> perspective(S left, S right, S bottom, S top, S zNear, S zFar,
                         bool directX = false)
{
    S l = left, r = right, b = bottom, t = top, n = zNear, f = zFar,
    rl = r - l, tb = t - b, nf = n - f;
    bool d = directX;
    return Mat<S, 4>((2 * n) / rl, 0, (r + l) / rl, 0,
                     0, (2 * n) / tb, (t + b) / tb, 0,
                     0, 0, ((d ? 0 : n) + f) / nf, ((d ? 1 : 2) * n * f) / nf,
                     0, 0, -1, 0);
}

/// Constructs a perspective view matrix.
template <class S>
IL Mat<S, 4> perspective(S fovYRadians, S aspect, S zNear, S zFar,
                         bool directX = false)
{
    S n = zNear, f = zFar, ys = 1 / tanf(0.5 * fovYRadians), nf = n - f;
    bool d = directX;
    return Mat<S, 4>(ys / aspect, 0, 0, 0,
                     0, ys, 0, 0,
                     0, 0, ((d ? 0 : n) + f) / nf, ((d ? 1 : 2) * n * f) / nf,
                     0, 0, -1, 0);
}

/// Constructs a perspective view matrix.
template <class S>
IL Mat<S, 4> infinitePerspective(S left, S right, S bottom, S top, S zNear,
                                 bool directX = false)
{
    S rl = right - left, tb = top - bottom;
    return Mat<S, 4>((2 * zNear) / rl, 0, (right + left) / rl, 0,
                     0, (2 * zNear) / tb, (top + bottom) / tb, 0,
                     0, 0, -1, (directX ? -zNear : -2 * zNear),
                     0, 0, -1, 0);
}

#undef min
#undef max

/// Returns a matrix with the component wise minimums
template <class S, U R, U C,
class D0, class D1, bool IR0, bool IR1, bool IV0, bool IV1>
IL Mat<S, R, C> min(const tMatBase<S, R, C, D0, IR0, IV0> &n,
                    const tMatBase<S, R, C, D1, IR1, IV1> &m)
{ Mat<S, R, C> r; F(j, C) F(i, R) r(i,j) = n(i,j) < m(i,j) ? n(i,j) : m(i,j);
    return r; }

/// Returns a matrix with the component wise minimums
template <class S, U R, U C, class D, bool IR, bool IV, class T>
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, Mat<S, R, C> >::type
min(const tMatBase<S, R, C, D, IR, IV> &m, const T &t)
{ Mat<S, R, C> r; F(j, C) F(i, R) r(i,j) = m(i,j) < t ? m(i,j) : t; 	return r; }

/// Returns a matrix with the component wise minimums
template <class S, U R, U C, class D, bool IR, bool IV, class T>
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, Mat<S, R, C> >::type
min(const T &t, const tMatBase<S, R, C, D, IR, IV> &m)
{ Mat<S, R, C> r; F(j, C) F(i, R) r(i,j) = t < m(i,j) ? t : m(i,j); 	return r; }

/// Returns a matrix with the component wise maximums
template <class S, U R, U C,
class D0, class D1, bool IR0, bool IR1, bool IV0, bool IV1>
IL Mat<S, R, C> max(const tMatBase<S, R, C, D0, IR0, IV0> &n,
                    const tMatBase<S, R, C, D1, IR1, IV1> &m)
{ Mat<S, R, C> r; F(j, C) F(i, R) r(i,j) = n(i,j) > m(i,j) ? n(i,j) : m(i,j);
    return r; }

/// Returns a matrix with the component wise maximums
template <class S, U R, U C, class D, bool IR, bool IV, class T>
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, Mat<S, R, C> >::type
max(const tMatBase<S, R, C, D, IR, IV> &m, const T &t)
{ Mat<S, R, C> r; F(j, C) F(i, R) r(i,j) = m(i,j) > t ? m(i,j) : t; return r; }

/// Returns a matrix with the component wise maximums
template <class S, U R, U C, class D, bool IR, bool IV, class T>
IL typename tMatUtil::EnableIf<tMatUtil::IsAri<T>::value, Mat<S, R, C> >::type
max(const T &t, const tMatBase<S, R, C, D, IR, IV> &m)
{ Mat<S, R, C> r; F(j, C) F(i, R) r(i,j) = t > m(i,j) ? t : m(i,j); return r; }

// End -------------------------------------------------------------------------

#undef DEF_MAT_BLOCK_E_WISE_OP_EQ_TILDA
#undef DEF_MAT_BLOCK_E_WISE_OP_EQ_SCALAR
#undef DEF_MAT_BLOCK_E_WISE_OP_EQ
#undef DEF_MAT_BLOCK_COMMON
#undef DEF_VEC_COMMON
#undef DEF_SQ_MAT_COMMON
#undef DEF_SQ_MAT_BLOCK_COMMON
#undef DEF_VEC

#undef LG
#undef DEF_MI
#undef VI
#undef MI
#undef F
#undef IL
#undef U

#define DEF_NAMESPACE };
DEF_NAMESPACE
#undef DEF_NAMESPACE

#pragma pop_macro("DI")
#pragma pop_macro("DEF_SQ_MAT_BLOCK_COMMON")
#pragma pop_macro("DEF_MAT_SINGLE_OP")
#pragma pop_macro("DEF_SIMPLE_MAT_MAT_OP")
#pragma pop_macro("DEF_MAT_SINGLE_IDENTITY")
#pragma pop_macro("LG")
#pragma pop_macro("DEF_MAT_BLOCK_COMMON")
#pragma pop_macro("DEF_MAT_SINGLE")
#pragma pop_macro("DEF_MI")
#pragma pop_macro("DEF_VEC")
#pragma pop_macro("DEF_VEC_SINGLE")
#pragma pop_macro("F")
#pragma pop_macro("VI")
#pragma pop_macro("DEF_MAT_BLOCK_E_WISE_OP_EQ_TILDA")
#pragma pop_macro("DEF_VEC_COMMON")
#pragma pop_macro("DEF_VEC_BLOCK")
#pragma pop_macro("DEF_NAMESPACE")
#pragma pop_macro("U")
#pragma pop_macro("IL")
#pragma pop_macro("DEF_SSE_INV")
#pragma pop_macro("SX")
#pragma pop_macro("DEF_MAT_BLOCK_E_WISE_OP_EQ")
#pragma pop_macro("SS")
#pragma pop_macro("DEF_MAT_TILDA_OP")
#pragma pop_macro("MI")
#pragma pop_macro("ST")
#pragma pop_macro("DEF_MAT_SCALAR_OP")
#pragma pop_macro("DEF_IS_ARI")
#pragma pop_macro("SI")
#pragma pop_macro("SM")
#pragma pop_macro("DEF_SQ_MAT_COMMON")
#pragma pop_macro("DEF_MAT_BLOCK_E_WISE_OP_EQ_SCALAR")
#pragma pop_macro("DEF_MAT_SINGLE_OP_SCALAR")

#endif
