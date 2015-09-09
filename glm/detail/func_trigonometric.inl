///////////////////////////////////////////////////////////////////////////////////
/// OpenGL Mathematics (glm.g-truc.net)
///
/// Copyright (c) 2005 - 2015 G-Truc Creation (www.g-truc.net)
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
/// 
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
/// 
/// Restrictions:
///		By making use of the Software for military purposes, you choose to make
///		a Bunny unhappy.
/// 
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
/// THE SOFTWARE.
///
/// @ref core
/// @file glm/detail/func_trigonometric.inl
/// @date 2008-08-03 / 2015-09-08
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#include "_vectorize.hpp"
#include <cmath>
#include <limits>
#include "../common.hpp"
#include "../exponential.hpp"
#include "../gtc/constants.hpp"

namespace glm
{
	// radians
	template <typename genType>
	GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType radians(genType degrees)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'radians' only accept floating-point input");

		return degrees * static_cast<genType>(0.01745329251994329576923690768489);
	}

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER GLM_CONSTEXPR vecType<T, P> radians(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(radians, v);
	}
	
	// degrees
	template <typename genType>
	GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType degrees(genType radians)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'degrees' only accept floating-point input");

		return radians * static_cast<genType>(57.295779513082320876798154814105);
	}

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER GLM_CONSTEXPR vecType<T, P> degrees(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(degrees, v);
	}


	// cos
	template <typename genType>
	GLM_FUNC_QUALIFIER genType cos_52s(genType const & x)
	{
		genType const xx(x * x);
		return (static_cast<genType>(0.9999932946) + xx * (static_cast<genType>(-0.4999124376) + xx * (static_cast<genType>(0.0414877472) + xx * static_cast<genType>(-0.0012712095))));
	}

#	if !GLM_FORCE_FLOAT_DETERMINISM
		using std::cos;
#	else
		template <typename genType>
		GLM_FUNC_QUALIFIER genType cos(genType x)
		{
			GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'cos' only accept floating-point input");

			auto wrap_angle = [](genType angle)
			{
				genType const before_two_pi(nextafterf(two_pi<genType>(), genType(0)));
				return abs<genType>(mod<genType>(angle, before_two_pi));
			};

			genType const result(wrap_angle(x));
			if (result<half_pi<genType>()) return cos_52s(result);
			if (result<pi<genType>()) return -cos_52s(pi<genType>() - result);
			if (result<(genType(3) * half_pi<genType>())) return -cos_52s(result - pi<genType>());
			return cos_52s(two_pi<genType>() - result);
		}
#	endif

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> cos(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(cos, v);
	}

	// sin
#	if !GLM_FORCE_FLOAT_DETERMINISM
		using ::std::sin;
#	else
		template <typename genType>
		GLM_FUNC_QUALIFIER genType sin(genType x)
		{
			return cos<genType>(half_pi<genType>() - x);
		}
#	endif

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> sin(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(sin, v);
	}

	// tan
	template <typename genType>
	GLM_FUNC_QUALIFIER genType tan_56s(genType const & x)
	{
		genType const xx(x * x);
		return (x*(static_cast<genType>(-3.16783027) + static_cast<genType>(0.134516124) * xx) / (static_cast<genType>(-4.033321984) + xx));
	}

#	if !GLM_FORCE_FLOAT_DETERMINISM
		using std::tan;
#	else
		template <typename genType>
		GLM_FUNC_QUALIFIER genType tan(genType x)
		{
			auto wrap_angle = [](genType angle)
			{
				genType const before_two_pi(nextafterf(two_pi<genType>(), genType(0)));
				return abs<genType>(mod<genType>(angle, before_two_pi));
			};

			genType result(wrap_angle(x));
			int const octant(static_cast<int>(result/quarter_pi<genType>()));
			assert(0 <= octant);
			assert(octant <= 7);
			if(0 == octant)
				return tan_56s(result * four_over_pi<genType>());
			else if(1 == octant)
				return genType(1) / tan_56s((half_pi<genType>() - result) * four_over_pi<genType>());
			else if(2 == octant)
				return genType(-1) / tan_56s((result - half_pi<genType>()) * four_over_pi<genType>());
			else if(3 == octant)
				return -tan_56s((pi<genType>() - result) * four_over_pi<genType>());
			else if(4 == octant)
				return tan_56s((result - pi<genType>()) * four_over_pi<genType>());
			else if(5 == octant)
				return genType(1) / tan_56s((three_over_two_pi<genType>() - result) * four_over_pi<genType>());
			else if(6 == octant)
				return genType(-1) / tan_56s((result - three_over_two_pi<genType>()) * four_over_pi<genType>());
			else if(7 == octant)
				return -tan_56s((two_pi<genType>() - result) * four_over_pi<genType>());
			else
				return x;
		}
#	endif

	// atan
	template <typename genType>
	GLM_FUNC_QUALIFIER genType atan(genType const & y, genType const & x)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'atan' only accept floating-point input");

#		ifdef GLM_FORCE_FLOAT_DETERMINISM
			genType sgn = sign(y) * sign(x);
			return abs(atan(y / x)) * sgn;
#		else
			return ::std::atan2(y, x);
#		endif
	}

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> atan(vecType<T, P> const & a, vecType<T, P> const & b)
	{
		return detail::functor2<T, P, vecType>::call(atan2, a, b);
	}

	template <typename genType>
	GLM_FUNC_QUALIFIER genType atan_66s(genType const & x)
	{
		genType const xx(x * x);
		return (x * (static_cast<genType>(1.6867629106) + static_cast<genType>(0.4378497304) * xx) / (static_cast<genType>(1.6867633134) + xx));
	}

#	if !GLM_FORCE_FLOAT_DETERMINISM
		using std::atan;
#	else
		template <typename genType>
		GLM_FUNC_QUALIFIER genType atan(genType x)
		{
			genType result = x;
			bool complement = false;
			bool region = false;
			bool sign = false;
			genType sixth_pi = (genType)0.52359877559829887307710723054658381403286156656252;
			genType tan_twelfth_pi = (genType)0.26794919243112270647255365849412763305719474618962;
			genType tan_sixth_pi = (genType)0.57735026918962576450914878050195745564760175127013;
			if (result < genType(0)) {
				result = -result;
				sign = true;
			}
			if (result > genType(1)){
				result = genType(1) / result;
				complement = true;
			}
			if (result > tan_twelfth_pi) {
				result = (result - tan_sixth_pi) / (genType(1) + tan_sixth_pi * result);
				region = true;
			}
			result = atan_66s(result);
			if (region) result += sixth_pi;
			if (complement) result = half_pi<genType>() - result;
			if (sign) result = -result;
			return result;
		}
#	endif

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> atan(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(atan, v);
	}

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> tan(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(tan, v);
	}

	// asin
#	if !GLM_FORCE_FLOAT_DETERMINISM
		using std::asin;
#	else
		template <typename genType>
		GLM_FUNC_QUALIFIER genType asin(genType x)
		{
			return atan(x / sqrt(static_cast<genType>(1) - x * x));
		}
#	endif

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> asin(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(asin, v);
	}

	// acos
#	if !GLM_FORCE_FLOAT_DETERMINISM
		using std::acos;
#	else
		template <typename genType>
		GLM_FUNC_QUALIFIER genType acos(genType x)
		{
			return half_pi<genType>() - asin(x);
		}
#	endif

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> acos(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(acos, v);
	}

	// sinh
	using std::sinh;

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> sinh(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(sinh, v);
	}

	// cosh
	using std::cosh;

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> cosh(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(cosh, v);
	}

	// tanh
	using std::tanh;

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> tanh(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(tanh, v);
	}

	// asinh
#	if GLM_HAS_CXX11_STL && !GLM_FORCE_FLOAT_DETERMINISM
		using std::asinh;
#	else
		template <typename genType> 
		GLM_FUNC_QUALIFIER genType asinh(genType const & x)
		{
			GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'asinh' only accept floating-point input");

			return (x < static_cast<genType>(0) ? static_cast<genType>(-1) : (x > static_cast<genType>(0) ? static_cast<genType>(1) : static_cast<genType>(0))) * log(abs(x) + sqrt(static_cast<genType>(1) + x * x));
		}
#	endif

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> asinh(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(asinh, v);
	}

	// acosh
#	if GLM_HAS_CXX11_STL && !GLM_FORCE_FLOAT_DETERMINISM
		using std::acosh;
#	else
		template <typename genType> 
		GLM_FUNC_QUALIFIER genType acosh(genType const & x)
		{
			GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'acosh' only accept floating-point input");

			if(x < static_cast<genType>(1))
				return static_cast<genType>(0);
			return log(x + sqrt(x * x - static_cast<genType>(1)));
		}
#	endif

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> acosh(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(acosh, v);
	}

	// atanh
#	if GLM_HAS_CXX11_STL && !GLM_FORCE_FLOAT_DETERMINISM
		using std::atanh;
#	else
		template <typename genType>
		GLM_FUNC_QUALIFIER genType atanh(genType const & x)
		{
			GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'atanh' only accept floating-point input");
		
			if(abs(x) >= static_cast<genType>(1))
				return 0;
			return static_cast<genType>(0.5) * log((static_cast<genType>(1) + x) / (static_cast<genType>(1) - x));
		}
#	endif

	template <typename T, precision P, template <typename, precision> class vecType>
	GLM_FUNC_QUALIFIER vecType<T, P> atanh(vecType<T, P> const & v)
	{
		return detail::functor1<T, T, P, vecType>::call(atanh, v);
	}
}//namespace glm
