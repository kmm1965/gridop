#pragma once

#include "cyclic_index.hpp"

template<typename T, class OP>
struct forward_first_derivative_operator_x : cyclic_index_operator_x<OP>
{
	typedef cyclic_index_operator_x<OP> super;
	typedef T value_type;

protected:
    constexpr forward_first_derivative_operator_x(typename super::index_type const& index, value_type h) : super(index), h(h){}

	template<typename EOP>
    __DEVICE
    constexpr get_value_type_t<EOP> operator()(size_t i_, const EOP &eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = i.first, iy = i.second;
        return (eobj_proxy[super::m_index[typename super::index_value_type(ix + 1, iy)]] -
			eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]]) / h;
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
    constexpr backward_first_derivative_operator_x(typename super::index_type const& index, value_type h) : super(index), h(h){}

	template<typename EOP>
    __DEVICE
    constexpr get_value_type_t<EOP> operator()(size_t i_, const EOP &eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = i.first, iy = i.second;
        return (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]] -
			eobj_proxy[super::m_index[typename super::index_value_type(ix - 1, iy)]]) / h;
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
    constexpr central_first_derivative_operator_x(typename super::index_type const& index, value_type h) : super(index), h2(h * 2){}

	template<typename EOP>
    __DEVICE
    constexpr get_value_type_t<EOP> operator()(size_t i_, EOP &eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = i.first, iy = i.second;
        return (eobj_proxy[super::m_index[typename super::index_value_type(ix + 1, iy)]] -
			eobj_proxy[super::m_index[typename super::index_value_type(ix - 1, iy)]]) / h2;
    }

private:
	const value_type h2;
};

template<typename T, typename TAG_X>
struct first_derivative_operator_x;

template<typename T>
struct first_derivative_operator_x<T, tag_main> : forward_first_derivative_operator_x<T, first_derivative_operator_x<T, tag_main> >
{
	typedef forward_first_derivative_operator_x<T, first_derivative_operator_x> super;

    template<typename TAG>
    struct get_tag_type;

    template<typename tag_y>
    struct get_tag_type<std::tuple<tag_main, tag_y> >
    {
        typedef std::tuple<tag_aux, tag_y> type;
    };

    constexpr first_derivative_operator_x(typename super::index_type const& index, T h) : super(index, h) {};

    IMPLEMENT_MATH_EVAL_OPERATOR(first_derivative_operator_x)
    REIMPLEMENT_OPERATOR_EXEC()
};

template<typename T>
struct first_derivative_operator_x<T, tag_aux> : backward_first_derivative_operator_x<T, first_derivative_operator_x<T, tag_aux> >
{
	typedef backward_first_derivative_operator_x<T, first_derivative_operator_x> super;

    template<typename TAG>
    struct get_tag_type;

    template<typename tag_y>
    struct get_tag_type<std::tuple<tag_aux, tag_y> >
    {
        typedef std::tuple<tag_main, tag_y> type;
    };

    constexpr first_derivative_operator_x(typename super::index_type const& index, T h) : super(index, h){};

    IMPLEMENT_MATH_EVAL_OPERATOR(first_derivative_operator_x)
    REIMPLEMENT_OPERATOR_EXEC()
};

/////////////
template<typename T, class OP>
struct forward_first_derivative_operator_y : cyclic_index_operator_y<OP>
{
	typedef cyclic_index_operator_y<OP> super;
	typedef T value_type;

protected:
    constexpr forward_first_derivative_operator_y(typename super::index_type const& index, value_type h) : super(index), h(h) {}

    template<typename EOP>
    __DEVICE
    constexpr get_value_type_t<EOP> operator()(size_t i_, const EOP& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = i.first, iy = i.second;
        return (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy + 1)]] -
			eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]]) / h;
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
    constexpr backward_first_derivative_operator_y(typename super::index_type const& index, value_type h) : super(index), h(h) {}

	template<typename EOP>
    __DEVICE
    constexpr get_value_type_t<EOP> operator()(size_t i_, const EOP& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = i.first, iy = i.second;
        return (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]] -
			eobj_proxy[super::m_index[typename super::index_value_type(ix, iy - 1)]]) / h;
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
    constexpr central_first_derivative_operator_y(value_type h) : h2(h * 2) {}

	template<typename EOP>
    __DEVICE
    constexpr get_value_type_t<EOP> operator()(size_t i_, EOP& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = i.first, iy = i.second;
        return (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy + 1)]] -
			eobj_proxy[super::m_index[typename super::index_value_type(ix, iy - 1)]]) / h2;
    }

private:
	const value_type h2;
};

template<typename T, typename TAG_Y>
struct first_derivative_operator_y;

template<typename T>
struct first_derivative_operator_y<T, tag_main> : forward_first_derivative_operator_y<T, first_derivative_operator_y<T, tag_main> >
{
	typedef forward_first_derivative_operator_y<T, first_derivative_operator_y> super;

    template<typename TAG>
    struct get_tag_type;

    template<typename tag_x>
    struct get_tag_type<std::tuple<tag_x, tag_main> >
    {
        typedef std::tuple<tag_x, tag_aux> type;
    };

    constexpr first_derivative_operator_y(typename super::index_type const& index, T h) : super(index, h) {};

    IMPLEMENT_MATH_EVAL_OPERATOR(first_derivative_operator_y)
    REIMPLEMENT_OPERATOR_EXEC()
};

template<typename T>
struct first_derivative_operator_y<T, tag_aux> : backward_first_derivative_operator_y<T, first_derivative_operator_y<T, tag_aux> >
{
	typedef backward_first_derivative_operator_y<T, first_derivative_operator_y> super;

    template<typename TAG>
    struct get_tag_type;

    template<typename tag_x>
    struct get_tag_type<std::tuple<tag_x, tag_aux> >
    {
        typedef std::tuple<tag_x, tag_main> type;
    };

    constexpr first_derivative_operator_y(typename super::index_type const& index, T h) : super(index, h) {};

    IMPLEMENT_MATH_EVAL_OPERATOR(first_derivative_operator_y)
    REIMPLEMENT_OPERATOR_EXEC()
};
