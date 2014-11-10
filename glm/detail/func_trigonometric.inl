///////////////////////////////////////////////////////////////////////////////////
/// OpenGL Mathematics (glm.g-truc.net)
///
/// Copyright (c) 2005 - 2014 G-Truc Creation (www.g-truc.net)
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
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
/// THE SOFTWARE.
///
/// @ref core
/// @file glm/core/func_trigonometric.inl
/// @date 2008-08-03 / 2011-06-15
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#include "_vectorize.hpp"
#include <cmath>
#include <limits>

namespace glm
{
	// radians
	template <typename genType>
	GLM_FUNC_QUALIFIER genType radians
	(
		genType const & degrees
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'radians' only accept floating-point input");

		return degrees * genType(0.01745329251994329576923690768489);
	}

	VECTORIZE_VEC(radians)
	
	// degrees
	template <typename genType>
	GLM_FUNC_QUALIFIER genType degrees
	(
		genType const & radians
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'degrees' only accept floating-point input");

		return radians * genType(57.295779513082320876798154814105);
	}

	VECTORIZE_VEC(degrees)

	template <typename genType>
	GLM_FUNC_QUALIFIER genType wrapAngle(genType const & angle)
	{
		genType const before_two_pi(nextafterf(two_pi<genType>(), genType(0)));
		return abs<genType>(mod<genType>(angle, before_two_pi));
	}

	VECTORIZE_VEC(wrapAngle)

	template <typename genType>
	GLM_FUNC_QUALIFIER genType cos_52s(genType const &  x)
	{
		genType const xx(x * x);
		return (genType(0.9999932946) + xx * (genType(-0.4999124376) + xx * (genType(0.0414877472) + xx * genType(-0.0012712095))));
	}

	VECTORIZE_VEC(cos_52s)

	// cos
	template <typename genType>
	GLM_FUNC_QUALIFIER genType cos(genType const & angle)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'cos' only accept floating-point input");
#ifdef GLM_FORCE_FLOAT_DETERMINISM
		genType const result(wrapAngle<genType>(angle));
		if (result<half_pi<genType>()) return cos_52s(result);
		if (result<pi<genType>()) return -cos_52s(pi<genType>() - result);
		if (result<(genType(3) * half_pi<genType>())) return -cos_52s(angle - pi<genType>());
		return cos_52s(two_pi<genType>() - result);
#else
		return genType(::std::cos(angle));
#endif
	}

	VECTORIZE_VEC(cos)

	// sin
	template <typename genType>
	GLM_FUNC_QUALIFIER genType sin
	(
		genType const & angle
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'sin' only accept floating-point input");
#ifdef GLM_FORCE_FLOAT_DETERMINISM
		return cos<genType>(half_pi<genType>() - angle);
#else
		return genType(::std::sin(angle));
#endif
	}

	VECTORIZE_VEC(sin)

	template <typename genType>
	GLM_FUNC_QUALIFIER genType tan_56s(genType const &  x)
	{
		genType const xx(x * x);
		return (x*(genType(-3.16783027) + genType(0.134516124) * xx) / (genType(-4.033321984) + xx));
	}

	// tan
	template <typename genType>
	GLM_FUNC_QUALIFIER genType tan
	(
		genType const & angle
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'tan' only accept floating-point input");
#ifdef GLM_FORCE_FLOAT_DETERMINISM
		genType const result(wrapAngle<genType>(angle));
		int const octant(result/quarter_pi<genType>());
		switch (octant){
		case 0: return tan_56s(result * four_over_pi<genType>());
		case 1: return genType(1) / tan_56s((half_pi<genType>() - result) * four_over_pi<genType>());
		case 2: return genType(-1) / tan_56s((result - half_pi<genType>()) * four_over_pi<genType>());
		case 3: return -tan_56s((pi<genType>() - result) * four_over_pi<genType>());
		case 4: return tan_56s((result - pi<genType>()) * four_over_pi<genType>());
		case 5: return genType(1) / tan_56s((three_over_two_pi<genType>() - result) * four_over_pi<genType>());
		case 6: return genType(-1) / tan_56s((result - three_over_two_pi<genType>()) * four_over_pi<genType>());
		case 7: return -tan_56s((two_pi<genType>() - result) * four_over_pi<genType>());
		}
#else
		return genType(::std::tan(angle));
#endif
	}

	VECTORIZE_VEC(tan)

	// asin
	template <typename genType>
	GLM_FUNC_QUALIFIER genType asin
	(
		genType const & x
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'asin' only accept floating-point input");
#ifdef GLM_FORCE_FLOAT_DETERMINISM
		return atan(x / sqrt(genType(1) - x*x));
#else
		return genType(::std::asin(x));
#endif
	}

	VECTORIZE_VEC(asin)

	// acos
	template <typename genType>
	GLM_FUNC_QUALIFIER genType acos
	(
		genType const & x
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'acos' only accept floating-point input");
#ifdef GLM_FORCE_FLOAT_DETERMINISM
		return half_pi<genType>() - asin(x);
#else
		return genType(::std::acos(x));
#endif
	}

	VECTORIZE_VEC(acos)

	// atan
	template <typename genType>
	GLM_FUNC_QUALIFIER genType atan
	(
		genType const & y, 
		genType const & x
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'atan' only accept floating-point input");
#ifdef GLM_FORCE_FLOAT_DETERMINISM
		genType sgn = sign(y) * sign(x);
		return abs(atan(y / x)) * sgn;
#else
		return genType(::std::atan2(y, x));
#endif
	}

	VECTORIZE_VEC_VEC(atan)

	template <typename genType>
	GLM_FUNC_QUALIFIER genType atan_66s(genType const & x)
	{
		genType const xx(x * x);
		return (x*(genType(1.6867629106) + genType(0.4378497304) * xx) / (genType(1.6867633134) + xx));
	}

	template <typename genType>
	GLM_FUNC_QUALIFIER genType atan
	(
		genType const & x
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'atan' only accept floating-point input");
#ifdef GLM_FORCE_FLOAT_DETERMINISM
		genType result = x;
		bool complement = false;
		bool region = false;
		bool sign = false;
		genType twelfth_pi(0.26179938779914943653855361527329190701643078328126);
		genType sixth_pi(0.52359877559829887307710723054658381403286156656252);
		genType tan_twelfth_pi(0.26794919243112270647255365849412763305719474618962);
		genType tan_sixth_pi(0.57735026918962576450914878050195745564760175127013);
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
#else
		return genType(::std::atan(x));
#endif
	}

	VECTORIZE_VEC(atan)

	// sinh
	template <typename genType> 
	GLM_FUNC_QUALIFIER genType sinh
	(
		genType const & angle
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'sinh' only accept floating-point input");

		return genType(std::sinh(angle));
	}

	VECTORIZE_VEC(sinh)

	// cosh
	template <typename genType> 
	GLM_FUNC_QUALIFIER genType cosh
	(
		genType const & angle
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'cosh' only accept floating-point input");

		return genType(std::cosh(angle));
	}

	VECTORIZE_VEC(cosh)

	// tanh
	template <typename genType>
	GLM_FUNC_QUALIFIER genType tanh
	(
		genType const & angle
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'tanh' only accept floating-point input");

		return genType(std::tanh(angle));
	}

	VECTORIZE_VEC(tanh)

	// asinh
	template <typename genType> 
	GLM_FUNC_QUALIFIER genType asinh
	(
		genType const & x
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'asinh' only accept floating-point input");
		
		return (x < genType(0) ? genType(-1) : (x > genType(0) ? genType(1) : genType(0))) * log(abs(x) + sqrt(genType(1) + x * x));
	}

	VECTORIZE_VEC(asinh)

	// acosh
	template <typename genType> 
	GLM_FUNC_QUALIFIER genType acosh
	(
		genType const & x
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'acosh' only accept floating-point input");

		if(x < genType(1))
			return genType(0);
		return log(x + sqrt(x * x - genType(1)));
	}

	VECTORIZE_VEC(acosh)

	// atanh
	template <typename genType>
	GLM_FUNC_QUALIFIER genType atanh
	(
		genType const & x
	)
	{
		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'atanh' only accept floating-point input");
		
		if(abs(x) >= genType(1))
			return 0;
		return genType(0.5) * log((genType(1) + x) / (genType(1) - x));
	}

	VECTORIZE_VEC(atanh)

}//namespace glm
