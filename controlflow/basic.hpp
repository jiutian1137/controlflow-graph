#pragma once

/// Basic of Control Flow with Graph.
///@license Free 
///@review 2022-9-23 
///@author LongJiangnan, Jiang1998Nan@outlook.com 
#define _CONTROLFLOW_BASIC_
#define _CONTROLFLOW_NAMED_SLOT 

#include <cassert>
#include <utility>
#include <memory>
#ifndef _STD_ALIGN_SIZE
#define _STD_ALIGN_SIZE
namespace std {
	inline size_t align(size_t _Bound, const size_t _Size) {
		//assert( _Bound != 0 && (_Bound & (_Bound - 1)) == 0 );
		return (_Size + (_Bound - 1)) & (~(_Bound - 1));
	}
}
#endif

#include <vector>

//#include <iosfwd>
#include <iostream>

#include <string>// to source|binary"opencl,glsl,hlsl" 

namespace controlflow {
	using invoke_result_t = size_t;
	constexpr invoke_result_t
		invoke_push_current = 0,
		invoke_pop_current  = 1,
		invoke_pop_last     = 1<<1,
		invoke_throw        = 1<<2,
		invoke_next         = 1<<3,
		invoke_nextpos      = 1<<4,
		invoke_nextpos_shift = 5,
		invoke_complete = invoke_pop_current|invoke_next,
		invoke_switch   = invoke_pop_current|invoke_next|invoke_nextpos,
		invoke_if_true  = invoke_switch,
		invoke_if_false = invoke_switch|(1<<invoke_nextpos_shift),
		invoke_continue = invoke_pop_current,
		invoke_break    = invoke_pop_current|invoke_pop_last,
		invoke_for      = invoke_push_current|invoke_next;

	struct vertex_property : public std::enable_shared_from_this<vertex_property> {
		virtual ~vertex_property() {}

		///@note should invoke before each frame, example:forloop all vertices .
		virtual void clear() noexcept {}
		/// while (completed) {
			///@return can skip is "false".
			virtual bool setup() { clear(); return true; }
			///@return end of loop is "false".
			virtual invoke_result_t invoke() { return controlflow::invoke_complete; }
#if 0
			/// used for control-flow, some result_slots cannot through, another result_slots can through
			///@note we use invoke_result_t instead this function.
			virtual bool next(size_t i) { return true; }
#endif
#if 0
			/// used for control-flow loop, always in the container(std::stack<..>) until end of loop.
			///@note control-flow loop use result by *::invoke().
			/// because we ignored the stack-pop process, it costs a lot, and
			///  the same loop order can be achieved by saving the current step loop-variables, except for the stack-pop process.
			virtual bool pop() { return true; }
#endif
		/// }

		///@note usually invoke ones, example:any graph search .
		/// ignore above two case functions if you want to implement "strict checking",
		/// example:{ ThrowIfFailed( dynamic_cast<const for_each_vertex_n&>(tg).attributes[i].length == .. ) }.
		virtual size_t num_arguments() const { return 0; }
		virtual bool set_argument(size_t i, const void*) { abort(); }
		virtual size_t num_results() const { return 0; }
		virtual const void* get_result(size_t i) const { abort(); }
		virtual bool link_argument(size_t i, const vertex_property& tg, size_t j) { return set_argument(i, tg.get_result(j)); }

#ifdef _CONTROLFLOW_NAMED_SLOT
		/// std::flat_map<std::string/*Name*/, size_t/*Index*/> argument_slot_mappings;
		/// std::flat_map<std::string/*Name*/, size_t/*Index*/> result_slot_mappings;
		/// some Name[s] will be constexpr, use these interfaces avoid dynamic memory in this case.
		///@note Even if there is no name, it can still be edited in the UI normally. 
		virtual size_t argument_index(const std::string&) const { abort(); }
		virtual std::string argument_name(size_t i) const { abort(); }
		virtual size_t result_index(const std::string&) const { abort(); }
		virtual std::string result_name(size_t i) const { abort(); }
#endif
	};

	using edge_property = std::pair< std::vector<std::pair<size_t/*source_slot*/,size_t/*target_slot*/>>, bool/*invoke*/ >;


	struct const_buffer_view {
		const unsigned char* _Ptr;
		size_t length;
		const unsigned char* data() const {
			return _Ptr;
		}
	};

	struct buffer_view : public const_buffer_view {
		unsigned char* data() const {
			return const_cast<unsigned char*>(const_buffer_view::_Ptr);
		}
	};

	struct type_info {
		size_t hash_code;/// scalar_name | ( (number_of_components - 1)<<bit_ceil(number_of_scalar_names) )
		size_t alignment;
		size_t length;///always aligned alignment unless 
	};

	struct member_type_info : type_info {
		size_t offset;
	};

	struct buffer : public vertex_property {
		std::vector<member_type_info> type_infos;
		struct in_block {
			unsigned char* variables;
		} out;

		buffer() noexcept = default;

		virtual ~buffer() override {
			if (out.variables != nullptr) {
				assert(!type_infos.empty());
				::operator delete(out.variables, std::align_val_t{ std::max(type_infos[0].alignment,
					__STDCPP_DEFAULT_NEW_ALIGNMENT__) });
			}
		}

		void _Construct_not_resize_typeinfos(const std::vector<controlflow::type_info>& type_infos) {
			assert(this->type_infos.size() == type_infos.size());
			if (!type_infos.empty()) {
				size_t max_length = std::align(type_infos[0].alignment, type_infos[0].length);
				for (size_t i = 1; i != type_infos.size(); ++i) {
					if (type_infos[i].alignment > type_infos[i - 1].alignment)
						max_length += (type_infos[i].alignment - type_infos[i - 1].alignment);
					max_length += std::align(type_infos[i].alignment, type_infos[i].length);
				}
				this->out.variables = reinterpret_cast<unsigned char*>(::operator new(max_length, std::align_val_t{ std::max(type_infos[0].alignment,
					__STDCPP_DEFAULT_NEW_ALIGNMENT__) }));

				size_t out_variables_i = reinterpret_cast<size_t>(this->out.variables);
				for (size_t i = 0; i != type_infos.size(); ++i) {
					out_variables_i = std::align(type_infos[i].alignment, out_variables_i);
					this->type_infos[i] = { type_infos[i], (out_variables_i - reinterpret_cast<size_t>(this->out.variables)) };
					out_variables_i += std::align(type_infos[i].alignment, type_infos[i].length);
				}
			}
		}

		explicit buffer(const std::vector<controlflow::type_info>& type_infos) : type_infos(type_infos.size()), out{ nullptr }
#ifdef _CONTROLFLOW_NAMED_SLOT
			, names(type_infos.size())
#endif
		{
			_Construct_not_resize_typeinfos(type_infos);
		}

		explicit buffer(std::vector<controlflow::type_info>&& type_infos) : buffer(static_cast<const std::vector<controlflow::type_info>&>(type_infos)) {}

		template<typename _Ty>
		_Ty& at(size_t i) {
			assert(i < num_results());
			return reinterpret_cast<_Ty&>(out.variables[type_infos[i].offset]);
		}

		template<typename _Ty>
		const _Ty& at(size_t i) const {
			assert(i < num_results());
			return reinterpret_cast<_Ty&>(out.variables[type_infos[i].offset]);
		}

		void move_from(size_t) {}

		template<typename _Ty0, typename... _TyN>
		void move_from(size_t i, _Ty0&& x0, _TyN&&... xN) {
			at<_Ty0>(i) = std::move(x0);
			return move_from(i + 1, std::forward<_TyN&&>(xN)...);
		}

		void copy_from(size_t) {}

		template<typename _Ty0, typename... _TyN>
		void copy_from(size_t i, const _Ty0& x0, _TyN&&... xN) {
			at<_Ty0>(i) = x0;
			return move_from(i + 1, std::forward<_TyN&&>(xN)...);
		}

		template<typename _Ty0, typename... _TyN>
		explicit buffer(_Ty0&& x0, _TyN&&... xN) : buffer(std::vector<controlflow::type_info>{
			controlflow::type_info{ typeid(_Ty0).hash_code(), alignof(_Ty0), sizeof(_Ty0) },
				controlflow::type_info{ typeid(_TyN).hash_code(), alignof(_TyN), sizeof(_TyN) }... }) {
			this->move_from(0, std::move(x0), std::forward<_TyN&&>(xN)...);
		}

		template<typename _Ty0, typename... _TyN>
		explicit buffer(const _Ty0& x0, _TyN&&... xN) : buffer(std::vector<controlflow::type_info>{
			controlflow::type_info{ typeid(_Ty0).hash_code(), alignof(_Ty0), sizeof(_Ty0) },
				controlflow::type_info{ typeid(_TyN).hash_code(), alignof(_TyN), sizeof(_TyN) }... }) {
			this->copy_from(1, x0, std::forward<_TyN&&>(xN)...);
		}

		virtual size_t num_results() const override {
			return type_infos.size();
		}

		virtual const void* get_result(size_t i) const override {
			assert(i < num_results());
			return &out.variables[type_infos[i].offset];
		}

#ifdef _CONTROLFLOW_NAMED_SLOT
		std::vector<std::string> names;

#if 0
		explicit buffer(const std::vector<std::pair<std::string, controlflow::type_info>>& name_info_s)
			: type_infos(type_infos.size()), out{ nullptr }, names(names.size()) {
			std::vector<controlflow::type_info> type_infos(name_info_s.size());
			_Construct_not_resize_typeinfos(type_infos);
			for (size_t i = 0; i != name_info_s.size(); ++i)
				names[i] = name_info_s[i].first;
		}

		explicit buffer(std::vector<std::pair<std::string, controlflow::type_info>>&& name_info_s)
			: type_infos(type_infos.size()), out{ nullptr }, names(names.size()) {
			std::vector<controlflow::type_info> type_infos(name_info_s.size());
			_Construct_not_resize_typeinfos(type_infos);
			for (size_t i = 0; i != name_info_s.size(); ++i)
				names[i] = name_info_s[i].first;
		}
#endif

		template<typename... _TyN>
		explicit buffer(const std::vector<std::string>& names, _TyN&&... xN) : type_infos(names.size()), out{ nullptr }, names(names.size()) {
			assert(names.size() == sizeof...(_TyN));
			_Construct_not_resize_typeinfos({ controlflow::type_info{ typeid(_TyN).hash_code(), alignof(_TyN), sizeof(_TyN) }... });
			this->move_from(0, std::forward<_TyN&&>(xN)...);
			for (size_t i = 0; i != names.size(); ++i)
				this->names[i] = names[i];
		}

		template<typename... _TyN>
		explicit buffer(std::vector<std::string>&& names, _TyN&&... xN) : buffer(static_cast<const std::vector<std::string>&>(names), std::forward<_TyN&&>(xN)...) {}

		virtual size_t result_index(const std::string& name) const override {
			for (const std::string& the_name : names)
				if (the_name == name)
					return std::distance(&names[0], &the_name);
			return static_cast<size_t>(-1);
		}

		virtual std::string result_name(size_t i) const override {
			assert(i < num_results());
			return names[i];
		}
#endif
	};

#if 0
#include <iostream>
#include "controlflow/basic.hpp"

	struct vec4 { float data[4]; };

	int main() {
		controlflow::buffer A(1, 2, 3, 4, 5, 6, 5.0f);
		for (size_t i = 0; i != A.num_results(); ++i) {
			std::cout << *((const int*)A.get_result(i)) << std::endl;
		}

		controlflow::buffer B = controlflow::buffer(std::vector<std::string>{
			"A", "B", "C", "D", "G" },
			10, 20, 30, vec4{},/*50UL,*/ 70ULL);
		for (size_t i = 0; i != B.num_results(); ++i) {
			std::cout << B.result_name(i) << ':' << *((const int*)B.get_result(i)) << std::endl;
		}

		return 0;
	}
#endif

	struct case_ : public controlflow::vertex_property {
		struct in_block {
			const size_t* witch_case;
		} in;

		virtual size_t num_arguments() const override { return 1; }

		virtual bool set_argument(size_t i, const void* arg) override {
			assert( i < num_arguments() );
			in.witch_case = (const size_t*)arg;
			return true;
		}
	};

	///@diagram
	/// +-------------------+
	/// |       switch_     |
	/// +-------------------+                  +---------+
	/// > witch_case  case0 > ---------------->|  case_  | -----> do thing0
	/// |             case1 > --------------\  +---------+
	/// |             case2 > -----------\  |  +---------+
	/// |             ...   >            |  +->|  case_  | -----> do thing1
	/// |             ...   >            |     +---------+
	///                                  |     +---------+
	///                                  +---->|  case_  | -----> do thing2
	///                                        +---------+
	///@note switch_.get_result(N) can only be set to "case_".
	struct switch_ : public controlflow::vertex_property {
		struct in_block {
			const size_t* witch_case;
		} in;

		virtual invoke_result_t invoke() override {
			return invoke_switch|(*in.witch_case);
		}

		virtual size_t num_arguments() const override { return 1; }

		virtual bool set_argument(size_t i, const void* arg) override {
			assert( i == num_arguments() );
			in.witch_case = (const size_t*)arg;
			return true;
		}
	};

	///@diagram
	/// +-----------------+
	/// |        if_      |
	/// +-----------------+         +---------+
	/// > condition  true > ------->|  case_  | ------> do "if" thing
	/// |                 |         +---------+
	/// |                 |         +---------+        +-----------------+
	/// |            false> ------->|  case_  | ------>|       if_       |
	/// +-----------------+         +---------+        +-----------------+         +---------+
	///                                                > condition  true > ------->|  case_  | ------> do "else if" thing
	///                                                |                 |         +---------+
	///                                                |                 |         +---------+
	///                                                |            false> ------->|  case_  | ------> do "else" thing
	///                                                +-----------------+         +---------+
	///@note if_.get_result(N) can only be set to "case_".
	struct if_ : public controlflow::vertex_property {
		struct in_block {
			const bool* condition;
		} in;
		struct out_block {
			size_t cases[2] = {0,1};
		} out;

		virtual invoke_result_t invoke() override {
			return (*in.condition) ? invoke_if_true : invoke_if_false;
		}

		virtual size_t num_arguments() const override { return 1; }

		virtual bool set_argument(size_t i, const void* arg) override {
			assert( i < num_arguments() );
			in.condition = (const bool*)arg;
			return true;
		}

		virtual size_t num_results() const override { return 2; }

		virtual const void* get_result(size_t i) const override {
			assert( i < num_results() );
			return &out.cases[i];
		}
	};

	struct for_ : public controlflow::vertex_property {
		size_t begin;
		size_t end;
		size_t index;
		struct in_block {// prior input variables
			const size_t* begin;
			const size_t* end;
		} in;
		struct out_block {
			size_t index;
		} out;

		for_() noexcept = default;

		explicit for_(size_t begin) : begin(begin), end(0), index(0), in{nullptr,nullptr}, out{0} {}
		
		explicit for_(size_t begin, size_t end) : begin(begin), end(end), index(0), in{nullptr,nullptr}, out{0} {}

		virtual void clear() noexcept override {
			index = 0;
		}

		virtual bool setup() override {
			if (in.begin != nullptr) 
				begin = *in.begin;
			if (in.end != nullptr) 
				end = *in.end;
			if (begin < end) {
				index = begin;
				return true;
			} else {
				return false;
			}
		}

		virtual invoke_result_t invoke() override {
			out.index = index++;
			return index != end ? invoke_for : invoke_complete;
		}

		virtual size_t num_arguments() const override { return 2; }

		virtual bool set_argument(size_t i, const void* arg) override {
			assert( i < num_arguments() );
			if (i == 0) { in.begin = (const size_t*)arg; }
			else if (i == 1) { in.end = (const size_t*)arg; }
			//else { in.is_break = (const bool*)arg; }
			return true;
		}

		virtual size_t num_results() const override { return 1; }

		virtual const void* get_result(size_t i) const override { assert( i < num_results() ); return &out.index; }

#ifdef _CONTROLFLOW_NAMED_SLOT
		virtual size_t argument_index(const std::string& name) const override { 
			return name == "begin" ? 0
				: name == "end" ? 1
				: name == "is_break" ? 2
				: static_cast<size_t>(-1);
		}

		virtual std::string argument_name(size_t i) const override { 
			assert( i < num_arguments() );
			return i == 0 ? "begin"
				: i == 1 ? "end"
				: "is_break";
		}

		virtual size_t result_index(const std::string& name) const override { return name == "index" ? 0 : static_cast<size_t>(-1); }

		virtual std::string result_name(size_t i) const override { assert( i < num_results() ); return "index"; }
#endif
	};


	template<typename _Ty>
	struct print : public controlflow::vertex_property {
		std::ostream* output;
		char endchar;
		struct in_block {
			const _Ty* value;
		} in;

		explicit print(std::ostream* output, char endchar = '\n') : output(output), endchar(endchar), in{nullptr} {}

		virtual size_t num_arguments() const override { return 1; }

		virtual bool set_argument(size_t i, const void* arg) override {
			assert( i < num_arguments() );
			in.value = (const _Ty*)arg;
			return true;
		}

		virtual invoke_result_t invoke() override { 
			if (in.value != nullptr)
				(*output) << (*in.value);
			(*output) << endchar;
			return invoke_complete;
		}
	};

#define _controlflow_binary_operator(NAME, OP) \
	template<typename _Ty1, typename _Ty2 = _Ty1> \
	struct NAME : public controlflow::vertex_property { \
		struct in_block { \
			const _Ty1* left; \
			const _Ty2* right; \
		} in; \
		struct out_block { \
			decltype( std::declval<_Ty1>() OP std::declval<_Ty2>() ) result; \
		} out; \
\
		virtual size_t num_arguments() const override { return 2; } \
\
		virtual bool set_argument(size_t i, const void* arg) override { \
			assert( i < num_arguments() ); \
			if (i == 0) in.left = (const _Ty1*)arg; \
			else in.right = (const _Ty2*)arg; \
			return true; \
		} \
\
		virtual size_t num_results() const override { return 1; } \
\
		virtual const void* get_result(size_t i) const override { \
			assert( i < num_results() ); \
			return &out.result; \
		} \
\
		virtual controlflow::invoke_result_t invoke() override { \
			out.result = (*in.left) ##OP## (*in.right); \
			return controlflow::invoke_complete; \
		} \
\
		virtual size_t argument_index(const std::string& name) const override { \
			return name == "left" ? 0 : 1; \
		} \
		virtual std::string argument_name(size_t i) const override { \
			assert( i < num_arguments() ); \
			return i == 0 ? "left" : "right"; \
		} \
		virtual size_t result_index(const std::string& name) const override { \
			return name == "result" ? 0 : static_cast<size_t>(-1); \
		} \
		virtual std::string result_name(size_t i) const override { \
			assert( i < num_results() ); \
			return "result"; \
		} \
	};

	_controlflow_binary_operator(plus, +)
	_controlflow_binary_operator(minus, -)
	_controlflow_binary_operator(multiplies, *)
	_controlflow_binary_operator(divides, /)
	_controlflow_binary_operator(modulus, %)
	_controlflow_binary_operator(equal_to, ==)
	_controlflow_binary_operator(not_equal_to, !=)

#define _controlflow_binary_assign_operator(NAME, ASSIGNOP) \
	template<typename _Ty1, typename _Ty2 = _Ty1> \
	struct NAME : public controlflow::vertex_property { \
		struct in_block { \
			_Ty1* left; \
			const _Ty2* right; \
		} in; \
		struct out_block { \
			_Ty1* result; \
		} out; \
\
		virtual size_t num_arguments() const override { return 2; } \
\
		virtual bool set_argument(size_t i, const void* arg) override { \
			assert( i < num_arguments() ); \
			if (i == 0) { assert(arg != nullptr); in.left = (_Ty1*)arg; } \
			else in.right = (const _Ty2*)arg; \
			return true; \
		} \
\
		virtual size_t num_results() const override { return 1; } \
\
		virtual const void* get_result(size_t i) const override { \
			assert( i < num_results() ); \
			return &out.result; \
		} \
\
		virtual controlflow::invoke_result_t invoke() override { \
			assert( in.left != nullptr ); \
			assert( in.right != nullptr ); \
			(*in.left) ##ASSIGNOP## (*in.right); \
			return controlflow::invoke_complete; \
		} \
\
		virtual size_t argument_index(const std::string& name) const override { \
			return name == "left" ? 0 : 1; \
		} \
		virtual std::string argument_name(size_t i) const override { \
			assert( i < num_arguments() ); \
			return i == 0 ? "left" : "right"; \
		} \
		virtual size_t result_index(const std::string& name) const override { \
			return name == "result" ? 0 : static_cast<size_t>(-1); \
		} \
		virtual std::string result_name(size_t i) const override { \
			assert( i < num_results() ); \
			return "result"; \
		} \
	};

	_controlflow_binary_assign_operator(assign, =)
	_controlflow_binary_assign_operator(assign_plus, +=)
	_controlflow_binary_assign_operator(assign_minus, -=)
	_controlflow_binary_assign_operator(assign_multiplies, *=)
	_controlflow_binary_assign_operator(assign_divides, /=)
	_controlflow_binary_assign_operator(assign_modulus, %=)
	_controlflow_binary_assign_operator(assign_bit_and, &=)
	_controlflow_binary_assign_operator(assign_bit_or, |=)
}