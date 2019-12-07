#pragma once

#include "cyclic_index.hpp"

template<class OP>
struct forward_average_operator : cyclic_index_operator<OP>
{
	typedef cyclic_index_operator<OP> super;

protected:
    constexpr forward_average_operator(typename super::index_type const& index) : super(index) {}

    template<typename EOP>
    constexpr get_value_type_t<EOP> operator()(size_t i, EOP const& eobj_proxy) const {
		return 0.5 * (eobj_proxy[i] + eobj_proxy[super::m_index[super::m_index.value(i) + 1]]);
	}
};

template<class OP>
struct backward_average_operator : cyclic_index_operator<OP>
{
    typedef cyclic_index_operator<OP> super;

protected:
    constexpr backward_average_operator(typename super::index_type const& index) : super(index) {}

    template<typename EOP>
    constexpr get_value_type_t<EOP> operator()(size_t i, EOP const& eobj_proxy) const {
		return 0.5 * (eobj_proxy[super::m_index[super::m_index.value(i) - 1]] + eobj_proxy[i]);
	}
};

template<typename TAG>
struct average_operator;

template<>
struct average_operator<tag_main> : forward_average_operator<average_operator<tag_main> >
{
    typedef forward_average_operator<average_operator> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
		static_assert(std::is_same<EO_TAG, tag_main>::value, "Tag types should be the same");
	};

    constexpr average_operator(typename super::index_type const& index) : super(index) {}

    IMPLEMENT_MATH_EVAL_OPERATOR(average_operator)
    REIMPLEMENT_OPERATOR_EXEC()
};

template<>
struct average_operator<tag_aux> : backward_average_operator<average_operator<tag_aux> >
{
    typedef backward_average_operator<average_operator> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
		static_assert(std::is_same<EO_TAG, tag_aux>::value, "Tag types should be the same");
	};

    constexpr average_operator(typename super::index_type const& index) : super(index) {}

    IMPLEMENT_MATH_EVAL_OPERATOR(average_operator)
    REIMPLEMENT_OPERATOR_EXEC()
};
