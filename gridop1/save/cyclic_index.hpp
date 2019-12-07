#pragma once

struct cyclic_index : index_base<cyclic_index>
{
	typedef cyclic_index type;
	typedef index_base<type> super;
	typedef size_t size_type;
	typedef isize_t index_type;

	cyclic_index(size_t size) : super(size){}

	__DEVICE __HOST
	size_t operator[](index_type i) const
	{
		const index_type N = (index_type)super::size();
		i %= N;
		assert(i > -N && i < N);
		if (i < 0)
			i += N;
		assert(i >= 0 && i < N);
		return i;
	}
};
