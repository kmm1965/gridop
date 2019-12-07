#pragma once

#include "cyclic_index.hpp"

template<class OP>
struct forward_average_operator_x : cyclic_index_operator_x<OP>
{
	typedef cyclic_index_operator_x<OP> super;

protected:
    forward_average_operator_x(typename super::index_type const& index) : super(index) {}

    template<typename EOP>
    __DEVICE
    typename EOP::value_type operator()(size_t i_, EOP const& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = std::get<0>(i), iy = std::get<1>(i);
		return 0.5 * (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]] +
			eobj_proxy[super::m_index[typename super::index_value_type(ix + 1, iy)]]);
	}
};

template<class OP>
struct backward_average_operator_x : cyclic_index_operator_x<OP>
{
	typedef cyclic_index_operator_x<OP> super;

protected:
    backward_average_operator_x(typename super::index_type const& index) : super(index) {}

    template<typename EOP>
    __DEVICE
	typename EOP::value_type operator()(size_t i_, EOP const& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = std::get<0>(i), iy = std::get<1>(i);
        return 0.5 * (eobj_proxy[super::m_index[typename super::index_value_type(ix - 1, iy)]] +
			eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]]);
    }
};

template<typename TAG_X>
struct average_operator_x;

template<>
struct average_operator_x<tag_main> : forward_average_operator_x<average_operator_x<tag_main> >
{
    typedef average_operator_x type;
	typedef forward_average_operator_x<type> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
		static_assert(std::is_same<std::tuple_element_t<0, EO_TAG>, tag_main>::value, "X tag types should be tag_main");
	};

	average_operator_x(typename super::index_type const& index) : super(index) {}

    REIMPLEMENT_OPERATOR_EXEC()
    REIMPLEMENT_MATH_EVAL_OPERATOR()
};

template<>
struct average_operator_x<tag_aux> : backward_average_operator_x<average_operator_x<tag_aux> >
{
    typedef average_operator_x type;
	typedef backward_average_operator_x<type> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
		static_assert(std::is_same<std::tuple_element_t<0, EO_TAG>, tag_aux>::value, "X tag type should be tag_aux");
	};

	average_operator_x(typename super::index_type const& index) : super(index) {}

    REIMPLEMENT_OPERATOR_EXEC()
    REIMPLEMENT_MATH_EVAL_OPERATOR()
};

template<class OP>
struct forward_average_operator_y : cyclic_index_operator_y<OP>
{
	typedef cyclic_index_operator_y<OP> super;

protected:
    forward_average_operator_y(typename super::index_type const& index) : super(index) {}

    template<typename EOP>
	typename EOP::value_type operator()(size_t i_, const EOP& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = std::get<0>(i), iy = std::get<1>(i);
        return 0.5 * (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]] 
                    + eobj_proxy[super::m_index[typename super::index_value_type(ix, iy + 1)]]);
    }
};

template<class OP>
struct backward_average_operator_y : cyclic_index_operator_y<OP>
{
	typedef cyclic_index_operator_y<OP> super;

protected:
	backward_average_operator_y(typename super::index_type const& index) : super(index) {}

	template<typename EOP>
	typename EOP::value_type operator()(index_type i_, const EOP& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = std::get<0>(i), iy = std::get<1>(i);
        return 0.5 * (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy - 1)]] 
                    + eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]]);
    }
};

template<typename TAG_Y>
struct average_operator_y;

template<>
struct average_operator_y<tag_main> : forward_average_operator_y<average_operator_y<tag_main> >
{
    typedef average_operator_y type;
	typedef forward_average_operator_y<type> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
		static_assert(std::is_same<std::tuple_element_t<1, EO_TAG>, tag_main>::value, "Y tag type should be tag_main");
	};

	average_operator_y(typename super::index_type const& index) : super(index) {}

    REIMPLEMENT_OPERATOR_EXEC()
    REIMPLEMENT_MATH_EVAL_OPERATOR()
};

template<>
struct average_operator_y<tag_aux> : backward_average_operator_y<average_operator_y<tag_aux> >
{
    typedef average_operator_y type;
	typedef backward_average_operator_y<type> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
		static_assert(std::is_same<std::tuple_element_t<1, EO_TAG>, tag_aux>::value, "Y tag type should be tag_aux");
	};

	average_operator_y(typename super::index_type const& index) : super(index) {}

    REIMPLEMENT_OPERATOR_EXEC()
    REIMPLEMENT_MATH_EVAL_OPERATOR()
};
