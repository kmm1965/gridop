#pragma once

#include "cyclic.hpp"

struct cyclic_index : index_base<cyclic_index>
{
	typedef isize_t value_type;

    constexpr cyclic_index(size_t size) : m_size(size){}

    constexpr size_t size() const {
        return m_size;
    }

    constexpr size_t operator[](value_type i) const
	{
        const value_type N = (value_type)m_size;
		i %= N;
		assert(i > -N && i < N);
		if (i < 0)
			i += N;
		assert(i >= 0 && i < N);
		return i;
	}

    constexpr value_type value(size_t i) const
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
    typedef cyclic_index index_type;
    typedef get_value_type_t<index_type> index_value_type;

	template<typename EO_TAG>
	struct get_tag_type
	{
		typedef other_tag_t<EO_TAG> type;
	};

protected:
    constexpr cyclic_index_operator(index_type const& index) : m_index(index.get_proxy()) {}

public:
    typename index_type::proxy_type const& index() const {
        return m_index;
    }

protected:
    const typename index_type::proxy_type m_index;
};
