#pragma once

#include "cyclic_index.hpp"

template<class OP>
struct forward_average_operator_x : cyclic_index_operator_x<OP>
{
    typedef cyclic_index_operator_x<OP> super;

protected:
    constexpr forward_average_operator_x(typename super::index_type const& index) : super(index) {}

    template<typename EOP>
    __DEVICE
    constexpr get_value_type_t<EOP> operator()(size_t i_, EOP const& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = i.first, iy = i.second;
		return 0.5 * (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]] +
			eobj_proxy[super::m_index[typename super::index_value_type(ix + 1, iy)]]);
	}
};

template<class OP>
struct backward_average_operator_x : cyclic_index_operator_x<OP>
{
    typedef cyclic_index_operator_x<OP> super;

protected:
    constexpr backward_average_operator_x(typename super::index_type const& index) : super(index) {}

    template<typename EOP>
    __DEVICE
    constexpr get_value_type_t<EOP> operator()(size_t i_, EOP const& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = i.first, iy = i.second;
        return 0.5 * (eobj_proxy[super::m_index[typename super::index_value_type(ix - 1, iy)]] +
			eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]]);
    }
};

template<typename TAG_X>
struct average_operator_x;

template<>
struct average_operator_x<tag_main> : forward_average_operator_x<average_operator_x<tag_main> >
{
	typedef forward_average_operator_x<average_operator_x> super;

    template<typename TAG>
    struct get_tag_type;

    template<typename tag_y>
    struct get_tag_type<std::tuple<tag_main, tag_y> >
    {
        typedef std::tuple<tag_aux, tag_y> type;
    };

    constexpr average_operator_x(typename super::index_type const& index) : super(index) {}

    IMPLEMENT_MATH_EVAL_OPERATOR(average_operator_x)
    REIMPLEMENT_OPERATOR_EXEC()
};

template<>
struct average_operator_x<tag_aux> : backward_average_operator_x<average_operator_x<tag_aux> >
{
	typedef backward_average_operator_x<average_operator_x> super;

    template<typename TAG>
    struct get_tag_type;

    template<typename tag_y>
    struct get_tag_type<std::tuple<tag_aux, tag_y> >
    {
        typedef std::tuple<tag_main, tag_y> type;
    };
    
    constexpr average_operator_x(typename super::index_type const& index) : super(index) {}

    IMPLEMENT_MATH_EVAL_OPERATOR(average_operator_x)
    REIMPLEMENT_OPERATOR_EXEC()
};

template<class OP>
struct forward_average_operator_y : cyclic_index_operator_y<OP>
{
    typedef cyclic_index_operator_y<OP> super;

protected:
    constexpr forward_average_operator_y(typename super::index_type const& index) : super(index) {}

    template<typename EOP>
    __DEVICE
    constexpr get_value_type_t<EOP> operator()(size_t i_, const EOP& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = i.first, iy = i.second;
        return 0.5 * (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]] +
            eobj_proxy[super::m_index[typename super::index_value_type(ix, iy + 1)]]);
    }
};

template<class OP>
struct backward_average_operator_y : cyclic_index_operator_y<OP>
{
    typedef OP type;
    typedef cyclic_index_operator_y<OP> super;

protected:
    constexpr backward_average_operator_y(typename super::index_type const& index) : super(index) {}

	template<typename EOP>
    __DEVICE
    constexpr get_value_type_t<EOP> operator()(size_t i_, const EOP& eobj_proxy) const
    {
        const typename super::index_value_type i = super::m_index.value(i_);
        const typename super::cyclic_index_value_type ix = i.first, iy = i.second;
        return 0.5 * (eobj_proxy[super::m_index[typename super::index_value_type(ix, iy - 1)]] +
            eobj_proxy[super::m_index[typename super::index_value_type(ix, iy)]]);
    }
};

template<typename TAG_Y>
struct average_operator_y;

template<>
struct average_operator_y<tag_main> : forward_average_operator_y<average_operator_y<tag_main> >
{
	typedef forward_average_operator_y<average_operator_y> super;

    template<typename TAG>
    struct get_tag_type;

    template<typename tag_x>
    struct get_tag_type<std::tuple<tag_x, tag_main> >
    {
        typedef std::tuple<tag_x, tag_aux> type;
    };

    constexpr average_operator_y(typename super::index_type const& index) : super(index) {}

    IMPLEMENT_MATH_EVAL_OPERATOR(average_operator_y)
    REIMPLEMENT_OPERATOR_EXEC()
};

template<>
struct average_operator_y<tag_aux> : backward_average_operator_y<average_operator_y<tag_aux> >
{
	typedef backward_average_operator_y<average_operator_y> super;

    template<typename TAG>
    struct get_tag_type;

    template<typename tag_x>
    struct get_tag_type<std::tuple<tag_x, tag_aux> >
    {
        typedef std::tuple<tag_x, tag_main> type;
    };

    constexpr average_operator_y(typename super::index_type const& index) : super(index) {}

    IMPLEMENT_MATH_EVAL_OPERATOR(average_operator_y)
    REIMPLEMENT_OPERATOR_EXEC()
};
