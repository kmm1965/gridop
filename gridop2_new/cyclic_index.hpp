#pragma once

#include "cyclic.hpp"

struct cyclic_index : index_base<cyclic_index>
{
    typedef isize_t value_type;

    cyclic_index(size_t size) : m_size(size){}

    __DEVICE __HOST
    size_t size() const {
        return m_size;
    }

    __DEVICE
    size_t operator[](value_type i) const
    {
        const value_type N = (value_type)m_size;
        i %= N;
        assert(i > -N && i < N);
        if (i < 0)
            i += N;
        assert(i >= 0 && i < N);
        return i;
    }

    __DEVICE
    value_type value(size_t i) const
    {
        assert(i < m_size);
        return i;
    }

private:
    const size_t m_size;
};

template<class MO>
struct cyclic_index_operator : math_operator<MO>
{
    typedef dim2_index<cyclic_index, cyclic_index> index_type;
    typedef typename index_type::value_type index_value_type;
    typedef typename cyclic_index::value_type cyclic_index_value_type;

protected:
    cyclic_index_operator(index_type const& index) : m_index(index.get_proxy()) {}

public:
    __DEVICE
    typename index_type::proxy_type const& index() const {
        return m_index;
    }

protected:
    const typename index_type::proxy_type m_index;
};

template<class MO>
struct cyclic_index_operator_x : cyclic_index_operator<MO>
{
    typedef cyclic_index_operator<MO> super;

	template<typename TAG>
	struct get_tag_type;

	template<typename tag_x, typename tag_y>
	struct get_tag_type<std::tuple<tag_x, tag_y> >
	{
		typedef std::tuple<other_tag_t<tag_x>, tag_y> type;
	};

protected:
    cyclic_index_operator_x(typename super::index_type const& index) : super(index) {}
};

template<class MO>
struct cyclic_index_operator_y : cyclic_index_operator<MO>
{
    typedef cyclic_index_operator<MO> super;

	template<typename TAG>
	struct get_tag_type;

	template<typename tag_x, typename tag_y>
	struct get_tag_type<std::tuple<tag_x, tag_y> >
	{
		typedef std::tuple<tag_x, other_tag_t<tag_y> > type;
	};

protected:
    cyclic_index_operator_y(typename super::index_type const& index) : super(index) {}
};
