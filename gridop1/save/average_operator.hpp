#pragma once

#include "cyclic.h"

template<typename TAG>
struct forward_average_operator : math_operator<TAG, forward_average_operator<TAG> >
{
	typedef forward_average_operator type;
	typedef math_operator<TAG, type> super;
	typedef isize_t index_type;

	template<typename EO_TAG>
	struct check_tag_type {
		static_assert(std::is_same<other_tag_t<TAG>, EO_TAG>::value, "Tags should be the same");
	};

	template<typename EOP>
	typename EOP::value_type operator()(index_type i, const EOP &eobj_proxy) const {
		return 0.5 * (eobj_proxy[i] + eobj_proxy[i + 1]);
	}

	REIMPLEMENT_MATH_EVAL_OPERATOR()
};

template<typename TAG>
struct backward_average_operator : math_operator<TAG, backward_average_operator<TAG> >
{
	typedef backward_average_operator type;
	typedef math_operator<TAG, type> super;
	typedef isize_t index_type;

	template<typename EO_TAG>
	struct check_tag_type {
		static_assert(std::is_same<other_tag_t<TAG>, EO_TAG>::value, "Tags should be the same");
	};

	template<typename EOP>
	typename EOP::value_type operator()(index_type i, const EOP &eobj_proxy) const {
		return 0.5 * (eobj_proxy[i - 1] + eobj_proxy[i]);
	}

	REIMPLEMENT_MATH_EVAL_OPERATOR()
};

template<typename TAG>
struct average_operator;

template<>
struct average_operator<tag_edge> : backward_average_operator<tag_edge> {};

template<>
struct average_operator<tag_cell> : forward_average_operator<tag_cell> {};

