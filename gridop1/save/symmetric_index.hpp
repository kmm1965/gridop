#pragma once

struct symmetric_index : index_base<symmetric_index>
{
	typedef symmetric_index type;
	typedef index_base<type> super;
	typedef size_t size_type;
	typedef isize_t index_type;

	symmetric_index(size_t size) : super(size){}

	__DEVICE __HOST
	size_t operator[](index_type i) const
	{
		const index_type N = (index_type)super::size();
		assert(i > -N && i < 2 * N - 1);
		if (i < 0) i = -i;
		else if (i >= N)
			i = N - (i - N + 2);
		assert(i >= 0 && i < N);
		return i;
	}
};
