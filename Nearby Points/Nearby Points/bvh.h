#ifndef bvh_335F0ABA_EA18_11E5_9CE9_5E5517507C66
#define bvh_335F0ABA_EA18_11E5_9CE9_5E5517507C66

// Here are the steps to use the BVH for accelerating
// searches for points within a certain radius.
//
// 1. Declare the BVH.
// For example, if you want a BVH for points,
// where each point is a 3D vector of double,
// you would do:
//
//     BVH<double, 3, true, PointId> bvh;
//
// Here, PointId can be something like an integer,
// or a custom class like:
//
//     struct PointId
//     {
//         size_t rodIndex;
//         size_t collocationPointIndex;
//     };
//
// 2. Add all the points. You can add them in a loop:
//
//     bvh.add(pointId, point);
//
// 3. Then, build the bvh:
//
//     bvh.build()
//
// 4. To get a points in a certain radius:
//
//     std::vector<PointId> intersectedPoints;
//     bvh.findPointsWithinRadius(intersectedPoints,
//                                queryPointCenter,
//                                queryPointRadius);
//
// 5. Then you can loop over the points and do stuff with them.
//
//     for (size_t i = 0; i < intersectedPoints.size(); ++i) {
//         PointId p = intersectedPoints[i];
//         // Do something with p...
//     }
//
// 6. If you want to clear the BVH,
// so that you can re-add entries and build it again,
// you can clear it with:
//
//     bvh.clear();
//

#include <vector>
#include <iostream>

#if defined(__SSE2__) || defined(_M_IX86_FP)
#include <emmintrin.h>
#endif

#pragma push_macro("F")
#pragma push_macro("align")
#pragma push_macro("DEF_FRIENDS")

#undef align
#undef F
#undef DEF_FRIENDS

#define F(I, UNTIL) for (unsigned I = 0; I < UNTIL; ++I)

#if defined(_MSC_VER)
#define align(x) __declspec(align(x))
#else
#if defined(__GNUC__)
#define align(x) __attribute__ ((aligned(x)))
#endif
#endif

namespace ActiveRodMesh
{
	/// BVH specialization for AABB / MBR.
	template <class Derived, bool ForPoints, class GeomId>
	class tBVHBase
	{
	public:
		
		/// Adds a geomId and its corresponding lower and upper bounds to the BVH.
		///
		/// The lower and upper bounds can be any vector objects that
		/// have the [] operator.
		///
		/// @param geomId  The geomId of the entry.
		/// @param lower   The lower bound of the AABB / MBR of the entru.
		/// @param upper   The upper bound of the AABB / MBR of the entry.
		template <class Vec>
		inline Derived & add(const GeomId &geomId, const Vec &lower, const Vec &upper)
		{ return static_cast<Derived*>(this)-> _addMBR(geomId, lower, upper); }
		
		/// Print out the points stored in the BVH. Useful for debugging.
		inline void printMBRs() const
		{ static_cast<const Derived *>(this)->_printMBRs(); }
	};
	
	/// BVH specialization for points.
	template <class Derived, class GeomId>
	class tBVHBase <Derived, true, GeomId>
	{
	public:
		
		/// Adds a geomId and its corresponding point to the BVH.
		///
		/// The point can be any vector object that has the [] operator.
		///
		/// @param geomId  The geomId of the entry.
		/// @param point   The coordinate of the entry's point.
		template <class Vec>
		inline Derived & add(const GeomId &geomId, const Vec &point)
		{ return static_cast<Derived *>(this)->_addPoint(geomId, point); }
		
		/// Find points that fall within a specified radius about a point.
		///
		/// The GeomIds of the entries which fall within the radius
		/// will be appended to the geomIds vector.
		///
		/// The point can be any vector object that has the [] operator.
		///
		/// The radius can be any integer / floating point value that can be
		/// converted to the Real type for this BVH.
		///
		/// It is recommended to use this function for such queries
		/// to make use of optimized SSE intrinsic tests.
		///
		/// @param geomIds  A std::vector to store the GeomIds of entries that
		///                 lie within the radius
		/// @param point    The point which the radius is about.
		/// @param radius   The radius.
		template <class Vec, class R, class GeomIdAllocator>
		inline void findPointsWithinRadius(std::vector<GeomId, GeomIdAllocator> &geomIds,
										   const Vec &point, const R &radius) const
		{
			return static_cast<const Derived *>(this)->
			_findPointsWithinRadius(geomIds, point, radius);
		}
		
		/// Print out the points stored in the BVH. Useful for debugging.
		inline void printPoints() const
		{ static_cast<const Derived *>(this)->_printPoints(); }
	};
	
	/// A Bounding Volume Hierarchy of Axis Aligned Bounding Boxes (AABB) /
	/// Minimum Bounding Rectangles (MBR)
	///
	/// @tparam Real      The type of Integer / Floating point value
	/// @tparam Dims      The number of dimensions of the space
	/// @tparam ForPoints Whether the BVH is to be specialized for points
	/// @tparam GeomId    A class to hold the index and data of entries.
	///                   Defaults to unsigned.
	template <class Real, unsigned Dims, bool ForPoints, class GeomId = unsigned>
	class BVH:
	public tBVHBase<BVH<Real, Dims, ForPoints, GeomId>, ForPoints, GeomId>
	{
		template <class, bool, class> friend class tBVHBase;
		
		template <class A> struct IsAri
		{
			template <int _, class A_> struct S { enum { v = 0 }; };
			template <int _> struct S <_, bool> { enum { v = 1 }; };
			template <int _> struct S <_, char> { enum { v = 1 }; };
			template <int _> struct S <_, unsigned char> { enum { v = 1 }; };
			template <int _> struct S <_, short> { enum { v = 1 }; };
			template <int _> struct S <_, unsigned short> { enum { v = 1 }; };
			template <int _> struct S <_, int> { enum { v = 1 }; };
			template <int _> struct S <_, unsigned int> { enum { v = 1 }; };
			template <int _> struct S <_, long> { enum { v = 1 }; };
			template <int _> struct S <_, unsigned long> { enum { v = 1 }; };
			template <int _> struct S <_, long long> { enum { v = 1 }; };
			template <int _> struct S <_, unsigned long long> { enum { v = 1 }; };
			template <int _> struct S <_, float> { enum { v = 1 }; };
			template <int _> struct S <_, double> { enum { v = 1 }; };
			template <int _> struct S <_, long double> { enum { v = 1 }; };
			enum { value = S <0, A>::v };
		};
		
		template <bool C, class S = void> struct EnableIf { };
		template <class S> struct EnableIf <true, S> { typedef S type; };
		
		template <int _, bool ChildInVec> class NodeBase;
		template <bool GeomIdInVec, bool IsPoint, bool UseSSE> class EntryBase;
		template <int _, class Real_, unsigned Dims_> struct MBRMBRIntersectionTest;
		
#define DEF_FRIENDS \
friend class BVH; \
template <int, bool> friend class NodeBase; \
template <bool, bool, bool> friend class EntryBase; \
template <int, class, unsigned> friend class MBRMBRIntersectionTest;
		
		template <bool UseSSE> class Vec
		{
			DEF_FRIENDS;
			
			template <class Real_, class V>
			static inline typename EnableIf<!IsAri<V>::value>::type
			_load(Real_ *d, const V &v)
			{ F(i, Dims) d[i] = v[i]; }
			
			template <class Real_, class V>
			static inline typename EnableIf<IsAri<V>::value>::type
			_load(Real_ *d, const V &v)
			{ F(i, Dims) d[i] = v; }
			
			template <int _, class Real_, unsigned Dims_, bool UseSSE_> struct Data
			{
				DEF_FRIENDS; friend class Vec;
				
				Real_ d[Dims_];
				
				inline void minInplace(const Data &o)
				{ F(i, Dims) d[i] = d[i] < o.d[i] ? d[i] : o.d[i]; }
				
				inline void maxInplace(const Data &o)
				{ F(i, Dims) d[i] = d[i] > o.d[i] ? d[i] : o.d[i]; }
				
				inline Real distSq(const Data &o) const
				{
					Real dSq = 0;
					F(i, Dims) dSq += (d[i] - o.d[i]) * (d[i] - o.d[i]);
					return dSq;
				}
			};
			
#if defined(__SSE2__) || defined(_M_IX86_FP)
			
			template <int _> class Data <_, float, 3, true>
			{
				DEF_FRIENDS; friend class Vec;
				
				float d[4];
				
				inline Data() {}
				
				inline Data(const Data &o) { _mm_storeu_ps(d, _mm_loadu_ps(o.d)); }
				
				inline Data & operator = (const Data &o)
				{ _mm_storeu_ps(d, _mm_loadu_ps(o.d)); return *this; }
				
				inline void minInplace(const Data &o)
				{
					__m128 m = _mm_min_ps(_mm_loadu_ps(d), _mm_loadu_ps(o.d));
					_mm_storeu_ps(d, m);
				}
				
				inline void maxInplace(const Data &o)
				{
					__m128 m = _mm_max_ps(_mm_loadu_ps(d), _mm_loadu_ps(o.d));
					_mm_storeu_ps(d, m);
				}
				
				inline float distSq(const Data &o) const
				{
					__m128 m = _mm_sub_ps(_mm_loadu_ps(d), _mm_loadu_ps(o.d));
					m = _mm_mul_ps(m, m);
					align(16) float p[4];
					_mm_storeu_ps(p, m);
					return p[0] + p[1] + p[2];
				}
			};
			
			template <int _> class Data <_, double, 3, true>
			{
				DEF_FRIENDS; friend class Vec;
				
				double d[4];
				
				inline Data() {}
				
				inline Data(const Data &o)
				{
					_mm_storeu_pd(d, _mm_loadu_pd(o.d));
					_mm_storeu_pd(d+2, _mm_loadu_pd(o.d+2));
				}
				
				inline Data & operator = (const Data &o)
				{
					_mm_storeu_pd(d, _mm_loadu_pd(o.d));
					_mm_storeu_pd(d+2, _mm_loadu_pd(o.d+2));
					return *this;
				}
				
				inline void minInplace(const Data &o)
				{
					__m128d m0 = _mm_min_pd(_mm_loadu_pd(d), _mm_loadu_pd(o.d));
					_mm_storeu_pd(d, m0);
					d[2] = d[2] < o.d[2] ? d[2] : o.d[2];
				}
				
				inline void maxInplace(const Data &o)
				{
					__m128d m0 = _mm_max_pd(_mm_loadu_pd(d), _mm_loadu_pd(o.d));
					_mm_storeu_pd(d, m0);
					d[2] = d[2] > o.d[2] ? d[2] : o.d[2];
				}
				
				inline double distSq(const Data &o) const
				{
					__m128d m0 = _mm_sub_pd(_mm_loadu_pd(d), _mm_loadu_pd(o.d));
					double m1 = d[2] - o.d[2];
					m0 = _mm_mul_pd(m0, m0);
					m1 = m1 * m1;
					align(16) double p[2];
					_mm_storeu_pd(p, m0);
					return p[0] + p[1] + m1;
				}
			};
#endif
			Data<0, Real, Dims, UseSSE> d;
			
		public:
			
			inline Real & operator[] (unsigned i) { return d.d[i]; }
			
			inline Real operator[] (unsigned i) const { return d.d[i]; }
			
			inline Vec() {}
			
			template <class V> inline Vec(const V &v) { _load(d.d, v); }
			
			template <class V> inline Vec & operator = (const V &v)
			{ _load(d.d, v); return *this; }
			
			inline Vec(const Vec &v): d(v.d) {}
			
			inline Vec & operator = (const Vec &b) { d = b.d; return *this; }
			
			template <class V> inline operator V () const volatile
			{ V v; F(i, Dims) v[i] = d.d[i]; return v; }
			
			inline void minInplace(const Vec &b) { d.minInplace(b.d); }
			
			inline void maxInplace(const Vec &b) { d.maxInplace(b.d); }
			
			inline Vec operator - (const Vec &o) const
			{ Vec v; F(i, Dims) v[i] = d.d[i] - o[i]; return v; }
			
			inline Vec operator + (const Vec &o) const
			{ Vec v; F(i, Dims) v[i] = d.d[i] + o[i]; return v; }
			
			inline Real distSq(const Vec &o) const
			{ return d.distSq(o.d); }
		};
		
		template <int _, bool ChildInVec> class NodeBase
		{
			Vec<true> b[2]; unsigned c[2];
		public:
			inline Vec<true> & lower() { return b[0]; }
			inline const Vec<true> & lower() const { return b[0]; }
			inline Vec<true> & upper() { return b[1]; }
			inline const Vec<true> & upper() const { return b[1]; }
			inline unsigned & child(unsigned i) { return c[i]; }
			inline const unsigned & child(unsigned i) const { return c[i]; }
		};
		
#if defined(__SSE2__) || defined(_M_IX86_FP)
		
		template <int _> class NodeBase <_, true>
		{
			Vec<true> b[2];
		public:
			inline Vec<true> & lower() { return b[0]; }
			inline const Vec<true> & lower() const { return b[0]; }
			inline Vec<true> & upper() { return b[1]; }
			inline const Vec<true> & upper() const { return b[1]; }
			inline unsigned & child(unsigned i)
			{ return *reinterpret_cast<unsigned *>(b[i].d.d + 3); }
			inline const unsigned & child(unsigned i) const
			{ return *reinterpret_cast<unsigned *>(b[i].d.d + 3); }
		};
		
#endif
		template <bool GeomIdInVec, bool IsPoint, bool UseSSE>
		class EntryBase
		{
			Vec<UseSSE> b[2 - IsPoint];
			GeomId g;
		public:
			inline GeomId & geomId() { return g; }
			inline const GeomId & geomId() const { return g; }
			inline Vec<UseSSE> & lower() { return b[0]; }
			inline const Vec<UseSSE> & lower() const { return b[0]; }
			inline Vec<UseSSE> & upper() { return b[!IsPoint]; }
			inline const Vec<UseSSE> & upper() const { return b[!IsPoint]; }
		};
		
#if defined(__SSE2__) || defined(_M_IX86_FP)
		
		template <bool IsPoint>
		class EntryBase <true, IsPoint, true>
		{
			Vec<true> b[2 - IsPoint];
		public:
			inline GeomId & geomId()
			{ return *reinterpret_cast<GeomId *>(b[0].d.d + 3); }
			inline GeomId geomId() const
			{ return *reinterpret_cast<GeomId *>(b[0].d.d + 3); }
			inline Vec<true> & lower() { return b[0]; }
			inline const Vec<true> & lower() const { return b[0]; }
			inline Vec<true> & upper() { return b[!IsPoint]; }
			inline const Vec<true> & upper() const { return b[!IsPoint]; }
		};
		
#endif
		template <int _, class Real_, unsigned Dims_>
		struct MBRMBRIntersectionTest
		{
			template <class M0, class M1>
			inline bool operator () (const M0 &m0, const M1 &m1)
			{
				F(i, Dims) if (m0.upper()[i] < m1.lower()[i]) return 0;
				F(i, Dims) if (m1.upper()[i] < m0.lower()[i]) return 0;
				return 1;
			}
		};
		
#if defined(__SSE2__) || defined(_M_IX86_FP)
		
		template <int _>
		struct MBRMBRIntersectionTest <_, float, 3>
		{
			template <class M0, class M1>
			inline bool operator () (const M0 &m0, const M1 &m1)
			{
				__m128 m0l = _mm_loadu_ps(m0.lower().d.d);
				__m128 m0u = _mm_loadu_ps(m0.upper().d.d);
				__m128 m1l = _mm_loadu_ps(m1.lower().d.d);
				__m128 m1u = _mm_loadu_ps(m1.upper().d.d);
				__m128 c0 = _mm_cmpgt_ps(m1l, m0u);
				__m128 c1 = _mm_cmpgt_ps(m0l, m1u);
				return !(_mm_movemask_ps(_mm_or_ps(c0, c1)) & 7);
			}
		};
		
		template <int _>
		struct MBRMBRIntersectionTest <_, double, 3>
		{
			template <class M0, class M1>
			inline bool operator () (const M0 &m0, const M1 &m1)
			{
				__m128d m0l0 = _mm_loadu_pd(m0.lower().d.d);
				__m128d m0u0 = _mm_loadu_pd(m0.upper().d.d);
				__m128d m1l0 = _mm_loadu_pd(m1.lower().d.d);
				__m128d m1u0 = _mm_loadu_pd(m1.upper().d.d);
				__m128d c00 = _mm_cmpgt_pd(m1l0, m0u0);
				__m128d c10 = _mm_cmpgt_pd(m0l0, m1u0);
				return (!_mm_movemask_pd(_mm_or_pd(c00, c10)) &&
						!(m1.lower().d.d[2] >= m0.upper().d.d[2] ||
						  m0.lower().d.d[2] >= m1.upper().d.d[2]));
			}
		};
#endif
		
#undef DEF_FRIENDS
		
#ifdef GL_VERSION
		
		template <class M>
		static inline void _drawMBR(const M &m)
		{
			const Vec<true> &l = m.lower();
			const Vec<true> &u = m.upper();
			
			glVertex3f(u[0], u[1], u[2]);
			glVertex3f(l[0], u[1], u[2]);
			
			glVertex3f(l[0], u[1], u[2]);
			glVertex3f(l[0], u[1], l[2]);
			
			glVertex3f(l[0], u[1], l[2]);
			glVertex3f(u[0], u[1], l[2]);
			
			glVertex3f(u[0], u[1], l[2]);
			glVertex3f(u[0], u[1], u[2]);
			
			glVertex3f(u[0], u[1], u[2]);
			glVertex3f(u[0], l[1], u[2]);
			
			glVertex3f(l[0], u[1], u[2]);
			glVertex3f(l[0], l[1], u[2]);
			
			glVertex3f(l[0], u[1], l[2]);
			glVertex3f(l[0], l[1], l[2]);
			
			glVertex3f(u[0], u[1], l[2]);
			glVertex3f(u[0], l[1], l[2]);
			
			glVertex3f(u[0], l[1], u[2]);
			glVertex3f(l[0], l[1], u[2]);
			
			glVertex3f(l[0], l[1], u[2]);
			glVertex3f(l[0], l[1], l[2]);
			
			glVertex3f(l[0], l[1], l[2]);
			glVertex3f(u[0], l[1], l[2]);
			
			glVertex3f(u[0], l[1], l[2]);
			glVertex3f(u[0], l[1], u[2]);
		}
		
		inline void _draw(unsigned i) const
		{
			if (i & 1) {
				_drawMBR(_leafs[i >> 1]);
			} else {
				_drawMBR(_nodes[i >> 1]);
				_draw(_nodes[i >> 1].child(0));
				_draw(_nodes[i >> 1].child(1));
			}
		}
#endif
		inline void _printMBRs() const
		{
			F(i, _numLeafs) {
				std::cout << "lower: [ ";
				F(j, 3) std::cout << _leafs[i].lower()[j] << " ";
				std::cout << "] T\n" << "upper: [ ";
				F(j, 3) std::cout << _leafs[i].upper()[j] << " ";
				std::cout << "] T\n\n";
			}
		}
		
		inline void _printPoints() const
		{
			F(i, _numLeafs) {
				std::cout << "[ ";
				F(j, 3) std::cout << _leafs[i].lower()[j] << " ";
				std::cout << "] T\n";
			}
		}
		
		inline void _grow()
		{
			unsigned oldNumSlots = _numSlots;
			Leaf *oldLeafs = _leafs;
			_numSlots <<= 1;
			_leafs = new Leaf[_numSlots];
			F(i, oldNumSlots) _leafs[i] = oldLeafs[i];
			delete [] oldLeafs;
		}
		
		
		typedef align(16) NodeBase
		<0, (sizeof(Real) >= sizeof(unsigned) && Dims == 3)> Node;
		
		typedef align(16) EntryBase
		<(sizeof(GeomId) <= sizeof(Real) && Dims == 3), ForPoints, true> Leaf;
		
		unsigned _numLeafs, _numSlots, _numNodes;
		unsigned _rootIndex;
		
		Node *_nodes;
		Leaf *_leafs;
		
		template <class Vec>
		inline BVH& _addMBR(const GeomId &geomId, const Vec &lower, const Vec &upper)
		{
			if (_numLeafs + 1 == _numSlots) _grow();
			_leafs[_numLeafs].lower() = lower;
			_leafs[_numLeafs].upper() = upper;
			_leafs[_numLeafs].geomId() = geomId;
			++_numLeafs;
			return *this;
		}
		
		template <class Vec>
		inline BVH& _addPoint(const GeomId &geomId, const Vec &point)
		{
			if (_numLeafs + 1 == _numSlots) _grow();
			_leafs[_numLeafs].lower() = point;
			_leafs[_numLeafs].geomId() = geomId;
			++_numLeafs;
			return *this;
		}
		
		template <class T>
		static inline void _swap(T &a, T &b) { T t = a; a = b; b = t; }
		
		template <class T, class Comp>
		static inline unsigned _sort3(T *t0, T *t1, T *t2, Comp &comp)
		{
			unsigned swaps = 0;
			if (comp(*t2, *t1)) { _swap(*t1, *t2); ++swaps; }
			if (comp(*t2, *t0)) { _swap(*t0, *t2); ++swaps; }
			if (comp(*t1, *t0)) { _swap(*t0, *t1); ++swaps; }
			return swaps;
		}
		
		template <class T, class Comp>
		static inline void _quickSelect(T *start, T *k, T *end, Comp comp)
		{
			while (true) {
			restart:
				if (k == end) return;
				switch (end - start) {
					case 0: case 1: return;
					case 2: if (comp(*--end, *start)) _swap(*start, *end); return;
					case 3: T *m = start; _sort3(start, ++m, --end, comp); return;
				}
				T *m = start + (end - start) / 2, *lm1 = end;
				unsigned swaps = _sort3(start, m, --lm1, comp);
				T *i = start, *j = lm1;
				if (!comp(*i, *m)) while (true) {
					if (i == --j) {
						++i; j = end;
						if (!comp(*start, *--j)) for (;; ++i) {
							if (i == j) return;
							if (comp(*start, *i)) { _swap(*i, *j); ++swaps; ++i; break; }
						}
						if (i == j) return;
						while (true) {
							while (!comp(*start, *i)) ++i;
							while (comp(*start, *--j)) ;
							if (i >= j) break;
							_swap(*i, *j); ++swaps; ++i;
						}
						if (k < i) return;
						start = i; goto restart;
					}
					if (comp(*j, *m)) { _swap(*i, *j); ++swaps; break; }
				}
				if (++i < j) for (;; ++i) {
					while (comp(*i, *m)) ++i;
					while (!comp(*--j, *m)) ;
					if (i >= j) break;
					_swap(*i, *j); ++swaps;
					if (m == i) m = j;
				}
				if (i != m && comp(*m, *i)) {
					_swap(*i, *m); ++swaps;
				}
				if (k == i) return;
				if (swaps == 0) {
					if (k < i) {
						for (j = m = start; ++j != i; m = j)
							if (comp(*j, *m)) goto notSorted;
						return;
					} else {
						for (j = m = i; ++j != end; m = j)
							if (comp(*j, *m)) goto notSorted;
						return;
					}
				}
			notSorted:
				if (k < i) end = i;
				else start = ++i;
			}
		}
		
		struct QuickSelectCompare
		{
			const unsigned axis;
			inline QuickSelectCompare(unsigned axis): axis(axis) {}
			inline bool operator() (const Leaf &a, const Leaf &b)
			{ return a.lower()[axis] < b.lower()[axis]; }
		};
		
		inline void _build(unsigned &index, unsigned begin, unsigned end)
		{
			if (end - begin == 1) {
				index = (begin << 1) | 1;
				return;
			}
			index = _numNodes << 1;
			Node &n = _nodes[_numNodes++];
			unsigned mid = begin + ((end - begin) >> 1);
			
			n.lower() = _leafs[begin].lower();
			n.upper() = _leafs[begin].upper();
			for (unsigned i = begin + 1; i < end; ++i) {
				n.lower().minInplace(_leafs[i].lower());
				n.upper().maxInplace(_leafs[i].upper());
			}
			
			if (end - begin > 2) {
				Vec<true> diff = n.upper() - n.lower();
				unsigned splitAxis = 0;
				for (unsigned i = 1; i < Dims; ++i)
					splitAxis = diff[i] > diff[splitAxis] ? i : splitAxis;
				
				_quickSelect(_leafs + begin, _leafs + mid, _leafs + end,
							 QuickSelectCompare(splitAxis));
				
			}
			
			_build(n.child(0), begin, mid);
			_build(n.child(1), mid, end);
		}
		
		template <class IntersectionTest, class IntersectionGeom, class GeomIdAllocator>
		inline void _findIntersections(std::vector<GeomId, GeomIdAllocator> &geomIds,
									   const IntersectionGeom &intersectionGeom,
									   const unsigned &index,
									   IntersectionTest &intersectionTest)
		{
			unsigned i = index >> 1;
			if (index & 1) {
				if (intersectionTest(intersectionGeom, _leafs[i].lower(), _leafs[i].upper()))
					geomIds.push_back(_leafs[i].geomId);
			} else if (intersectionTest(intersectionGeom, _nodes[i].lower(), _nodes[i].upper())) {
				F(j, 2) _findIntersections(geomIds, intersectionGeom, _nodes[i].child(j),
										   intersectionTest);
			}
		}
		
		template <class IntersectionTest, class IntersectionGeom, class GeomIdAllocator>
		inline void _findIntersections(std::vector<GeomId, GeomIdAllocator> &geomIds,
									   const IntersectionGeom &intersectionGeom,
									   IntersectionTest &intersectionTest)
		{ _findIntersections(geomIds, intersectionGeom, _rootIndex, intersectionTest); }
		
		template <class GeomIdAllocator>
		inline void _findMBRIntersections(std::vector<GeomId, GeomIdAllocator> &geomIds,
										  const Node &n, const unsigned &index)
		{
			MBRMBRIntersectionTest<0, Real, Dims> intersects;
			unsigned i = index >> 1;
			if (index & 1) {
				if (intersects(n, _leafs[i])) geomIds.push_back(_leafs[i].geomId());
			} else if (intersects(n, _nodes[i])) {
				F(j, 2) _findMBRIntersections(geomIds, n, _nodes[i].child(j));
			}
		}
		
		template <class GeomIdAllocator>
		inline void _findPointsWithinRadius(std::vector<GeomId, GeomIdAllocator> &geomIds,
											const Node &n,
											const Vec<true> &p,
											const Real &rSq,
											const unsigned &index) const
		{
			unsigned i = index >> 1;
			if (index & 1) {
				if (_leafs[i].lower().distSq(p) < rSq)
					geomIds.push_back(_leafs[i].geomId());
			} else if (MBRMBRIntersectionTest<0, Real, Dims>()(n, _nodes[i])) {
				F(j, 2) _findPointsWithinRadius(geomIds, n, p, rSq, _nodes[i].child(j));
			}
		}
		
		
		template <class V, class R, class GeomIdAllocator>
		inline void _findPointsWithinRadius(std::vector<GeomId, GeomIdAllocator> &geomIds,
											const V &point, const R &radius) const
		{
			Real rSq = radius * radius;
			Vec<true> p(point);
			Node n; n.lower() = p; n.upper() = p;
			F(i, Dims) n.lower()[i] -= radius;
			F(i, Dims) n.upper()[i] += radius;
			_findPointsWithinRadius(geomIds, n, p, rSq, _rootIndex);
		}
		
		inline void _copy(const BVH &o)
		{
			_numLeafs = o._numLeafs;
			_numNodes = o._numNodes;
			_nodes = new Node[_numNodes];
			F(i, _numNodes) _nodes[i] = o._nodes[i];
			_leafs = new Leaf[_numLeafs];
			F(i, _numLeafs) _leafs[i] = o._leafs[i];
		}
		
		inline void _release() { delete [] _nodes; delete [] _leafs; }
		
	public:
		
		/// Construct a copy of the other BVH.
		inline BVH(const BVH &o) { copy(o); }
		
		/// Copies the other BVH.
		inline BVH & operator = (const BVH &o) { _release(); copy(o); return *this; }
		
		/// Constructor.
		///
		/// @param initialCapacity  An optional parameter to hint how much
		///                         memory to preallocate for the BVH.
		inline BVH(size_t initialCapacity = 128)
		{
			_numLeafs = 0;
			for (_numSlots = 1; _numSlots < initialCapacity; _numSlots <<= 1) ;
			_nodes = NULL;
			_leafs = new Leaf[_numSlots];
		}
		
		/// Destructor.
		inline ~BVH() { _release(); }
		
		/// Builds the BVH.
		/// Please call this function before using the BVH to find intersections.
		inline BVH& build()
		{
			delete [] _nodes;
			_numNodes = 0;
			_nodes = new Node[_numLeafs - 1];
			_build(_rootIndex, 0, _numLeafs);
			return *this;
		}
		
		/// Find intersections with a custom geometry.
		///
		/// The GeomIds of the entries which have their Axis-Aligned Bounding Box (AABB) /
		/// Minimum Bounding Rectangle (MBR) intersected with the current geometry
		/// will be appended to the geomIds vector.
		///
		/// IntersectionTest can be a function class that implements the following:
		///
		/// @code
		/// bool operator () (const IntersectionGeom &geom,
		///                   const Vec& lowerBound,
		///                   const Vec& upperBound)
		/// @endcode
		///
		/// Vec can be any vector class that has the [] operator.
		/// The lower and upper bound of an node's AABB / MBR
		/// will be automatically converted to Vec and passed into the function.
		///
		/// The function should return true if the geometry intersect the AABB/ MBR
		/// to help determine which branches of the tree to search.
		///
		/// For points, the lower and upper bounds will be equal for the entries.
		///
		/// @param geomIds     A std::vector to store the GeomIds of intersected entries.
		/// @param geom        The geometry to be intersected with.
		/// @param intersects  The function class to test for intersection between the
		///                    custom geometry and AABB / MBR of the nodes and entries.
		template <class IntersectionTest, class IntersectionGeom, class GeomIdAllocator>
		inline void findIntersections(std::vector<GeomId, GeomIdAllocator> &geomIds,
									  const IntersectionGeom &geom,
									  IntersectionTest intersects) const
		{ const_cast<BVH *>(this)->_findIntersections(geomIds, geom, intersects); }
		
		/// Find intersections with a Axis-Aligned Bounding Box (AABB) /
		/// Minimum Bounding Rectangle (MBR) defined by lower and upper bounds.
		///
		/// The GeomIds of the entries which have their AABB/ MBR
		/// intersected will be appended to the geomIds vector.
		///
		/// The lower and upper bounds can be any vector objects that have
		/// the [] operator.
		///
		/// It is recommended to use this method to test for AABB / MBR
		/// intersections to make use of optimized SSE intersection tests.
		///
		/// @param geomIds  A std::vector to store the GeomIds of intersected entries.
		/// @param lower    The lower bound of the AABB / MBR to intersect with.
		/// @param upper    The upper bound of the AABB / MBR to intersect with.
		template <class V0, class V1, class GeomIdAllocator>
		inline void findMBRIntersections(std::vector<GeomId, GeomIdAllocator> &geomIds,
										 const V0 &lower, const V1 &upper) const
		{
			Node n; n.lower() = lower; n.upper() = upper;
			const_cast<BVH *>(this)->_findMBRIntersections(geomIds, n, _rootIndex);
		}
		
		/// Returns the number of entries in the BVH.
		inline unsigned size() const { return _numLeafs; }
		
		/// Clears the BVH.
		inline BVH & clear()
		{
			_numLeafs = 0;
			delete [] _nodes;
			_nodes = NULL;
			return *this;
		}
		
#ifdef GL_VERSION
		/// Draws the BVH with OpenGL.
		/// Requires an OpenGL context.
		inline void draw() const
		{
			glDisable(GL_LIGHTING);
			glBegin(GL_LINES);
			_draw(_rootIndex);
			glEnd();
			glEnable(GL_LIGHTING);
		}
#endif
		
#undef align
#undef F
#pragma pop_macro("F")
#pragma pop_macro("align")
#pragma pop_macro("DEF_FRIENDS")
		
	};
	
};

#endif
