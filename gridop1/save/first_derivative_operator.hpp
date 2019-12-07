#pragma once

#include "cyclic.h"

template<typename TAG, typename T>
struct forward_first_derivative_operator : math_operator<TAG, forward_first_derivative_operator<TAG, T> >
{
	typedef forward_first_derivative_operator type;
	typedef math_operator<TAG, type> super;
	typedef T value_type;
	typedef isize_t index_type;

	template<typename EO_TAG>
	struct check_tag_type {
		static_assert(std::is_same<other_tag_t<TAG>, EO_TAG>::value, "Tags should be the same");
	};

	forward_first_derivative_operator(value_type h) : h(h){}

	template<typename EOP>
	typename EOP::value_type operator()(index_type i, const EOP &eobj_proxy) const {
		return (eobj_proxy[i + 1] - eobj_proxy[i]) / h;
	}

	REIMPLEMENT_MATH_EVAL_OPERATOR()

private:
	const value_type h;
};

template<typename TAG, typename T>
struct backward_first_derivative_operator : math_operator<TAG, backward_first_derivative_operator<TAG, T> >
{
	typedef backward_first_derivative_operator type;
	typedef math_operator<TAG, type> super;
	typedef T value_type;
	typedef isize_t index_type;

	template<typename EO_TAG>
	struct check_tag_type {
		static_assert(std::is_same<other_tag_t<TAG>, EO_TAG>::value, "Tags should be the same");
	};

	backward_first_derivative_operator(value_type h) : h(h){}

	template<typename EOP>
	typename EOP::value_type operator()(index_type i, const EOP &eobj_proxy) const {
		return (eobj_proxy[i] - eobj_proxy[i - 1]) / h;
	}

	REIMPLEMENT_MATH_EVAL_OPERATOR()

private:
	const value_type h;
};

template<typename TAG, typename T>
struct central_first_derivative_operator : math_operator<TAG, central_first_derivative_operator<TAG, T> >
{
	typedef central_first_derivative_operator type;
	typedef math_operator<TAG, type> super;
	typedef T value_type;
	typedef isize_t index_type;

	template<typename EO_TAG>
	struct check_tag_type {
		static_assert(std::is_same<other_tag_t<TAG>, EO_TAG>::value, "Tags should be the same");
	};

	central_first_derivative_operator(value_type h) : h2(h * 2){}

	template<typename EOP>
	typename EOP::value_type operator()(index_type i, EOP &eobj_proxy) const {
		return (eobj_proxy[i + 1] - eobj_proxy[i - 1]) / h2;
	}

	REIMPLEMENT_MATH_EVAL_OPERATOR()

private:
	const value_type h2;
};

template<typename TAG>
struct first_derivative_operator;

template<>
struct first_derivative_operator<tag_edge> : backward_first_derivative_operator <tag_edge, double> 
{
    first_derivative_operator<tag_edge>(double h): backward_first_derivative_operator <tag_edge, double>(h){};
};

template<>
struct first_derivative_operator<tag_cell> : forward_first_derivative_operator <tag_cell, double> {
    first_derivative_operator<tag_cell>(double h): forward_first_derivative_operator <tag_cell, double>(h){};
};
