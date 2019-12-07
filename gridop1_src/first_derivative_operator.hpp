#pragma once

#include "cyclic_index.hpp"

template<typename T, class OP>
struct forward_first_derivative_operator : cyclic_index_operator<OP>
{
	typedef cyclic_index_operator<OP> super;
	typedef T value_type;
    
protected:
	constexpr forward_first_derivative_operator(typename super::index_type const& index, value_type h) : super(index), h(h){}

public:
    template<typename EOP>
    constexpr get_value_type_t<EOP> operator()(size_t i, const EOP &eobj_proxy) const {
        return (eobj_proxy[super::m_index[super::m_index.value(i) + 1]] - eobj_proxy[i]) / h;
    }

private:
	const value_type h;
};

template<typename T, class OP>
struct backward_first_derivative_operator : cyclic_index_operator<OP>
{
    typedef cyclic_index_operator<OP> super;
	typedef T value_type;

protected:
    constexpr backward_first_derivative_operator(typename super::index_type const& index, value_type h) :  super(index), h(h){}

public:
    template<typename EOP>
    constexpr get_value_type_t<EOP> operator()(size_t i, const EOP &eobj_proxy) const {
        return (eobj_proxy[i] - eobj_proxy[super::m_index[super::m_index.value(i) - 1]]) / h;
    }

private:
    const value_type h;
};

template<typename T, class OP>
struct central_first_derivative_operator : cyclic_index_operator<OP>
{
    typedef cyclic_index_operator<OP> super;
	typedef T value_type;

protected:
    constexpr central_first_derivative_operator(typename super::index_type const& index, value_type h) : super(index), h2(h * 2){}

public:
    template<typename EOP>
    constexpr get_value_type_t<EOP> operator()(size_t i_, EOP &eobj_proxy) const
    {
        const get_value_type_t<typename super::index_type> i = super::m_index.value(i_);
		return (eobj_proxy[super::m_index[i + 1]] - eobj_proxy[super::m_index[i - 1]]) / h2;
	}

private:
    const value_type h2;
};

template<typename T, typename TAG>
struct first_derivative_operator;

template<typename T>
struct first_derivative_operator<T, tag_main> : forward_first_derivative_operator<T, first_derivative_operator<T, tag_main> >
{
	typedef forward_first_derivative_operator<T, first_derivative_operator> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
		static_assert(std::is_same<tag_main, EO_TAG>::value, "Tag types should be the same");
	};

    constexpr first_derivative_operator(typename super::index_type const& index, T h) : super(index, h) {};

    IMPLEMENT_MATH_EVAL_OPERATOR(first_derivative_operator)
    REIMPLEMENT_OPERATOR_EXEC()
};

template<typename T>
struct first_derivative_operator<T, tag_aux> : backward_first_derivative_operator<T, first_derivative_operator<T, tag_aux> >
{
	typedef backward_first_derivative_operator<T, first_derivative_operator> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
		static_assert(std::is_same<tag_aux, EO_TAG>::value, "Tag types should be the same");
	};

    constexpr first_derivative_operator(typename super::index_type const& index, T h) : super(index, h){};

    IMPLEMENT_MATH_EVAL_OPERATOR(first_derivative_operator)
    REIMPLEMENT_OPERATOR_EXEC()
};
