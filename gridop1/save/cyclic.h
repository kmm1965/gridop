#pragma once

struct tag_cell;
struct tag_edge;

template<typename T>
struct other_tag;

template<>
struct other_tag<tag_cell>
{
	typedef tag_edge type;
};

template<>
struct other_tag<tag_edge>
{
	typedef tag_cell type;
};

template<typename T>
using other_tag_t = typename other_tag<T>::type;
