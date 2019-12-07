#pragma once

#include "../gridop1_src/cyclic_index.hpp"

template<class OP>
struct first_derivative_operator : cyclic_index_operator<OP>
{
    typedef cyclic_index_operator<OP> super;

    template<typename T>
    struct get_value_type;

    template<class Unit, class Y>
    struct get_value_type<boost::units::quantity<Unit, Y> >
    {
        typedef typename boost::units::divide_typeof_helper<boost::units::quantity<Unit, Y>, length_type>::type type;
    };

protected:
    constexpr first_derivative_operator(typename super::index_type const& index) : super(index) {}
};

template<class OP>
struct forward_first_derivative_operator : first_derivative_operator<OP>
{
	typedef first_derivative_operator<OP> super;
    
protected:
    constexpr forward_first_derivative_operator(typename super::index_type const& index, const length_type &h) : super(index), h(h){}

public:
    template<typename EOP>
    constexpr typename super::template get_value_type<get_value_type_t<EOP> >::type operator()(size_t i, const EOP &eobj_proxy) const {
        return (eobj_proxy[super::m_index[super::m_index.value(i) + 1]] - eobj_proxy[i]) / h;
    }

private:
    const length_type h;
};

template<class OP>
struct backward_first_derivative_operator : first_derivative_operator<OP>
{
    typedef first_derivative_operator<OP> super;
	
protected:
    constexpr backward_first_derivative_operator(typename super::index_type const& index, const length_type &h) :  super(index), h(h){}

public:
    template<typename EOP>
    constexpr typename super::template get_value_type<get_value_type_t<EOP> >::type operator()(size_t i, const EOP &eobj_proxy) const {
        return (eobj_proxy[i] - eobj_proxy[super::m_index[super::m_index.value(i) - 1]]) / h;
    }

private:
    const length_type h;
};

template<class OP>
struct central_first_derivative_operator : first_derivative_operator<OP>
{
    typedef OP type;
    typedef first_derivative_operator<OP> super;

protected:
    constexpr central_first_derivative_operator(typename super::index_type const& index, const length_type &h) : super(index), h2(h * 2.){}

public:
    template<typename EOP>
    constexpr typename super::template get_value_type<get_value_type_t<EOP> >::type operator()(size_t i_, EOP &eobj_proxy) const
    {
        const get_value_type_t<typename super::index_type> i = super::m_index.value(i_);
		return (eobj_proxy[super::m_index[i + 1]] - eobj_proxy[super::m_index[i - 1]]) / h2;
	}

private:
    const length_type h2;
};

template<typename TAG>
struct first_derivative_operator;

template<>
struct first_derivative_operator<tag_main> : forward_first_derivative_operator<first_derivative_operator<tag_main> >
{
	typedef forward_first_derivative_operator<first_derivative_operator> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
		static_assert(std::is_same<tag_main, EO_TAG>::value, "Tag types should be the same");
	};

    constexpr first_derivative_operator(typename super::index_type const& index, const length_type &h) : super(index, h) {};

    IMPLEMENT_MATH_EVAL_OPERATOR(first_derivative_operator)
    REIMPLEMENT_OPERATOR_EXEC()
};

template<>
struct first_derivative_operator<tag_aux> : backward_first_derivative_operator<first_derivative_operator<tag_aux> >
{
	typedef backward_first_derivative_operator<first_derivative_operator> super;

	template<typename EO_TAG>
	struct get_tag_type : super::template get_tag_type<EO_TAG>
	{
		static_assert(std::is_same<tag_aux, EO_TAG>::value, "Tag types should be the same");
	};

    constexpr first_derivative_operator(typename super::index_type const& index, const length_type &h) : super(index, h){};

    IMPLEMENT_MATH_EVAL_OPERATOR(first_derivative_operator)
    REIMPLEMENT_OPERATOR_EXEC()
};
