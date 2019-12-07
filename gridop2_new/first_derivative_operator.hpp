#pragma once

#include "cyclic_index.hpp"

template<typename T, class OP>
struct forward_first_derivative_operator_x : cyclic_index_operator_x<OP>
{
	typedef cyclic_index_operator_x<OP> super;
	typedef T value_type;

protected:
	forward_first_derivative_operator_x(typename super::index_type const& index, value_type h) : super(index), h(h){}

	template<typename EOP>
    __DEVICE
    typename EOP::value_type operator()(size_t i_, const EOP &eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = std::get<0>(i), iy = std::get<1>(i);
        return (  eobj_proxy[super::m_index[typename super::index_value_type(ix + 1, iy)]]
			    - eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]]) / h;
    }

private:
	const value_type h;
};

template<typename T, class OP>
struct backward_first_derivative_operator_x : cyclic_index_operator_x<OP>
{
	typedef cyclic_index_operator_x<OP> super;
	typedef T value_type;

protected:
	backward_first_derivative_operator_x(typename super::index_type const& index, value_type h) : super(index), h(h){}

	template<typename EOP>
	typename EOP::value_type operator()(size_t i_, const EOP &eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = std::get<0>(i), iy = std::get<1>(i);
        return (  eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]] 
			    - eobj_proxy[super::m_index[typename super::index_value_type(ix - 1, iy)]]) / h;
    }

private:
	const value_type h;
};

template<typename T, class OP>
struct central_first_derivative_operator_x : cyclic_index_operator_x<OP>
{
	typedef cyclic_index_operator_x<OP> super;
	typedef T value_type;

protected:
	central_first_derivative_operator_x(typename super::index_type const& index, value_type h) : super(index), h2(h * 2){}

	template<typename EOP>
	typename EOP::value_type operator()(size_t i_, EOP &eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = std::get<0>(i), iy = std::get<1>(i);
        return (eobj_proxy[super::m_index[super::index_value_type(ix + 1, iy)]] -
			eobj_proxy[super::m_index[super::index_value_type(ix - 1, iy)]]) / h2;
    }

private:
	const value_type h2;
};

template<typename T, typename TAG_X>
struct first_derivative_operator_x;

template<typename T>
struct first_derivative_operator_x<T, tag_main> : forward_first_derivative_operator_x<T, first_derivative_operator_x<T, tag_main> >
{
    typedef first_derivative_operator_x type;
	typedef forward_first_derivative_operator_x<T, type> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
        static_assert(std::is_same<std::tuple_element_t<0, EO_TAG>, tag_main>::value, "X tag types should be tag_main");
    };

	first_derivative_operator_x(typename super::index_type const& index, T h) : super(index, h) {};

    REIMPLEMENT_OPERATOR_EXEC()
    REIMPLEMENT_MATH_EVAL_OPERATOR()
};

template<typename T>
struct first_derivative_operator_x<T, tag_aux> : backward_first_derivative_operator_x<T, first_derivative_operator_x<T, tag_aux> >
{
    typedef first_derivative_operator_x type;
	typedef backward_first_derivative_operator_x<T, type> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
        static_assert(std::is_same<std::tuple_element_t<0, EO_TAG>, tag_aux>::value, "X tag types should be tag_aux");
    };

	first_derivative_operator_x(typename super::index_type const& index, T h) : super(index, h){};

    REIMPLEMENT_OPERATOR_EXEC()
    REIMPLEMENT_MATH_EVAL_OPERATOR()
};

/////////////
template<typename T, class OP>
struct forward_first_derivative_operator_y : cyclic_index_operator_y<OP>
{
	typedef cyclic_index_operator_y<OP> super;
	typedef T value_type;

protected:
	forward_first_derivative_operator_y(typename super::index_type const& index, value_type h) : super(index), h(h) {}

	template<typename EOP>
	typename EOP::value_type operator()(size_t i_, const EOP& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = std::get<0>(i), iy = std::get<1>(i);
        return (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy + 1)]] 
			  - eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]]) / h;
    }

private:
	const value_type h;
};

template<typename T, class OP>
struct backward_first_derivative_operator_y : cyclic_index_operator_y<OP>
{
	typedef cyclic_index_operator_y<OP> super;
	typedef T value_type;

protected:
	backward_first_derivative_operator_y(typename super::index_type const& index, value_type h) : super(index), h(h) {}

	template<typename EOP>
	typename EOP::value_type operator()(size_t i_, const EOP& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = std::get<0>(i), iy = std::get<1>(i);
        return (  eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]] 
			    - eobj_proxy[super::m_index[typename super::index_value_type(ix, iy - 1)]]) / h;
    }

private:
	const value_type h;
};

template<typename T, class OP>
struct central_first_derivative_operator_y : cyclic_index_operator_y<OP>
{
	typedef cyclic_index_operator_y<OP> super;
	typedef T value_type;

protected:
	central_first_derivative_operator_y(value_type h) : h2(h * 2) {}

	template<typename EOP>
	typename EOP::value_type operator()(size_t i_, EOP& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = std::get<0>(i), iy = std::get<1>(i);
        return (eobj_proxy[super::m_index[super::index_value_type(ix, iy + 1)]] -
			eobj_proxy[super::m_index[super::index_value_type(ix, iy - 1)]]) / h2;
    }

private:
	const value_type h2;
};

template<typename T, typename TAG_Y>
struct first_derivative_operator_y;

template<typename T>
struct first_derivative_operator_y<T, tag_main> : forward_first_derivative_operator_y<T, first_derivative_operator_y<T, tag_main> >
{
    typedef first_derivative_operator_y type;
	typedef forward_first_derivative_operator_y<T, type> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
        static_assert(std::is_same<std::tuple_element_t<1, EO_TAG>, tag_main>::value, "Y tag type should be tag_main");
    };

	first_derivative_operator_y(typename super::index_type const& index, double h) : super(index, h) {};

    REIMPLEMENT_OPERATOR_EXEC()
    REIMPLEMENT_MATH_EVAL_OPERATOR()
};

template<typename T>
struct first_derivative_operator_y<T, tag_aux> : backward_first_derivative_operator_y<T, first_derivative_operator_y<T, tag_aux> >
{
    typedef first_derivative_operator_y type;
	typedef backward_first_derivative_operator_y<T, type> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
        static_assert(std::is_same<std::tuple_element_t<1, EO_TAG>, tag_aux>::value, "Y tag type should be tag_aux");
    };

	first_derivative_operator_y(typename super::index_type const& index, double h) : super(index, h) {};

    REIMPLEMENT_OPERATOR_EXEC()
    REIMPLEMENT_MATH_EVAL_OPERATOR()
};
