#pragma once

/// Graphics Pipeline of Control Flow with Graph.
///@license Free 
///@review 2022-9-23 
///@author LongJiangnan, Jiang1998Nan@outlook.com 
#define _CONTROLFLOW_GRAPHICS_

#include "basic.hpp"

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
#if !(defined _HAS_CXX23 && _HAS_CXX23 != 0)
#include "../geometry/graph.hpp"
#else
#include <flat_map>
#endif

///@diagram
/// +--------------------------------------------------------------------------------------------+
/// | rasterization graphics pipeline /                                                          |
/// +--------------------------------/                                                           |
/// |                                                                                            |
/// |                   ---> .... --->                                                           |
/// |                 /                \                                                         |
/// |  +-------------+                  +----------------+                                       |
/// |  | for_each_vertex_n | ----> .... ----> | primitive |                                       |
/// |  +-------------+                  +----------------+                                       |
/// |                 \                /         |                                               |
/// |                   ---> .... --->           |                                               |
/// |                                            |                                               |
/// |              /----------<---------<--------/                                               |
/// |              |                                                                             |
/// |              |              ---> ... --->                                                  |
/// |              |            /               \                                                |
/// |  +-----------------------+                 +------------------+                            |
/// |  | raster_primitive      | ----> ... ----> | write_image_node |                            |
/// |  +-----------------------+                 +------------------+                            |
/// |                           \               /                                                |
/// |                             ---> ... --->                                                  |
/// |                                                                                            |
/// +--------------------------------------------------------------------------------------------+
namespace controlflow { namespace graphics {
	enum class primitive_types {
		none          = 0,
		point_list    = 1,
		line_list     = 2,
		triangle_list = 3,
			
		setup_one_vertex = (1<<2),
		setup_two_vertices = (2<<2),
		line_stripe    = line_list|setup_one_vertex,
		triangle_strip = triangle_list|setup_two_vertices,
		triangle_fan   = triangle_list|setup_two_vertices|(1<<4),
	};

	constexpr size_t _Get_primitive_vertices(primitive_types type) {
		return size_t(type) & 3;
	}

	constexpr size_t _Get_primitive_setup_vertices(primitive_types type) {
		return (size_t(type) >> 2) & 3;
	}

	static constexpr const char* vertex_index    = "vertex_index";
	static constexpr const char* primitive_index = "primitive_index";
	static constexpr const char* NDCposition     = "NDCposition";

	template<typename... attribute_types>
	struct for_each_vertex_n {};

	template<>
	struct for_each_vertex_n<> : public controlflow::vertex_property {
		size_t vertices;
		size_t index;
		struct attribute_binding {
			std::pair<attribute_type, size_t/*offset*/> type;
			std::pair<const_buffer_view, size_t/*stride*/> stream;
		};
		std::flat_map<std::string, attribute_binding> attributes;
		struct out_block {
			unsigned char* attributes;
			size_t vertex_index;
		} out;
	
		for_each_vertex_n() noexcept = default;

		for_each_vertex_n(size_t vertices, const std::vector<std::string>& names, const std::vector<attribute_type>& types, const std::vector<std::pair<const_buffer_view,size_t/*stride*/>>& streams)
			: vertices(vertices), index(0), attributes(), out{ nullptr, 0 } {
			assert( names.size() == types.size() && names.size() == streams.size() );
			if (!names.empty()) {
				size_t max_length = std::align(types[0].alignment, types[0].length);
				for (size_t i = 1; i != names.size(); ++i) {
					if (types[i].alignment > types[i - 1].alignment) 
						max_length += (types[i].alignment - types[i-1].alignment);
					max_length += std::align(types[i].alignment, types[i].length);
				}
				this->out.attributes = reinterpret_cast<unsigned char*>(::operator new(max_length, std::align_val_t{std::max(types[0].alignment,
					__STDCPP_DEFAULT_NEW_ALIGNMENT__)}));

				size_t out_attribute_i = reinterpret_cast<size_t>(this->out.attributes);
				for (size_t i = 0; i != names.size(); ++i) {
					out_attribute_i = std::align(types[i].alignment, out_attribute_i);
					this->attributes.insert_or_assign(names[i], 
						attribute_binding{ { types[i], reinterpret_cast<size_t>(this->out.attributes) - out_attribute_i }, streams[i] });
					out_attribute_i += std::align(types[i].alignment, types[i].length);
				}
			}
		}

		virtual ~for_each_vertex_n() override {
			if (out.attributes != nullptr) {
				assert( !attributes.empty() );
				::operator delete(out.attributes, std::align_val_t{std::max(attributes.values()[0].type.first.alignment,
					__STDCPP_DEFAULT_NEW_ALIGNMENT__)});
			}
		}

		virtual size_t num_results() const override {
			return attributes.size() + 1;
		}

		virtual size_t result_index(const std::string& name) const override {
			auto attribute_i = attributes.find(name);
			return attribute_i != attributes.end() ? std::distance(attributes.begin(), attribute_i)
				: "index" == name ? attributes.size()
				: static_cast<size_t>(-1);
		}

		virtual std::string result_name(size_t i) const override {
			assert( i < num_results() );
			return i < attributes.size() ? attributes.keys()[i]
				: "index";
		}

		virtual const void* get_result(size_t i) const override {
			assert( i < num_results() );
			return i < attributes.size() ? (const void*)std::next(out.attributes, attributes.values()[i].type.second)
				: (const void*)(&out.index);
		}

		virtual void clear() noexcept override {
			index = 0;
		}

		virtual bool setup() override {
			if (vertices == 0) return false;// skip null forloop.
			assert( index == vertices );// no exception when forloop completed.
			index = 0;
			return true;
		}

		virtual void invoke() override {
			for (auto first = attributes.values().begin(); first != attributes.values().end(); ++first) 
				std::copy_n(std::next(first->stream.first.data() + first->stream.second * index), first->type.first.length, 
					std::next(out.attributes, first->type.second));
			out.index = index++;
		}

		virtual bool pop() override {
			return index == vertices;
		}
	};

	template<typename... attribute_types>
	struct primitive : public controlflow::vertex_property {
		/// ...
	};

	template<>
	struct primitive<> : public controlflow::vertex_property {
		primitive_types type;
		size_t vertices;
		size_t setup_vertices;
		size_t index;
		std::pair<attribute_type,/*offset = 0,*/size_t/*stride*/>
			NDCposition;
		std::flat_map<std::string, std::pair<attribute_type,std::pair<size_t/*offset*/,size_t/*stride*/>>> 
			attributes;
		struct in_block {
			const unsigned char* 
				NDCposition;
			std::vector<const unsigned char*> 
				attributes;
		} in;
		struct out_block {
			/// [0]: | vertex[0].NDCposition, vertex[1].NDCposition, vertex[2].NDCposition    |
			/// [1]: | vertex[0].attributes0, vertex[1].attributes0, vertex[2].attributes0 |
			/// [2]: | vertex[0].attributes1, vertex[1].attributes1, vertex[2].attributes1 |
			/// [.]: | ...                    ...                    ...                   |
			/// [N]: | vertex[0].attributesN, vertex[1].attributesN, vertex[2].attributesN |
			///@note this matrix is transposed from original arrangement because faster cache (similar case is Math Matrix).
			unsigned char* 
				primitive_matrix;
		} out;

		primitive(primitive_types type, const std::vector<std::string>& names, const std::vector<attribute_type>& types, 
			const attribute_type& NDCposition_type = { 0, /*sizeof(vec4) = */16, /*alignof(vec4) = */16 })
			: type(type), vertices(_Get_primitive_vertices(type)), setup_vertices(_Get_primitive_setup_vertices(type)), index(0),
			NDCposition(NDCposition_type, 0), attributes(), in{ nullptr, std::vector<const unsigned char*>(names.size(),nullptr) }, out{ nullptr } {
			assert( names.size() == types.size() );
			size_t max_length = std::align(NDCposition_type.alignment, NDCposition_type.length) * vertices;
			if (!names.empty()) {
				for (size_t i = 0; i != names.size(); ++i) {
					if (types[i].alignment > types[i - 1].alignment) 
						max_length += (types[i].alignment - types[i-1].alignment);
					max_length += std::align(types[i].alignment, types[i].length) * vertices;
				}
			}
			this->out.primitive_matrix = reinterpret_cast<unsigned char*>(::operator new(max_length, std::align_val_t{std::max(NDCposition_type.alignment,
				__STDCPP_DEFAULT_NEW_ALIGNMENT__)}));

			NDCposition.second = std::align(NDCposition_type.alignment, NDCposition_type.length);
			if (!names.empty()) {
				size_t out_primitive_matrix_i = reinterpret_cast<size_t>(this->out.primitive_matrix) + NDCposition.second * vertices;
				for (size_t i = 0; i != names.size(); ++i) {
					out_primitive_matrix_i = std::align(types[i].alignment, out_primitive_matrix_i);
					this->attributes.insert_or_assign(names[i], 
						std::pair(types[i], std::pair(reinterpret_cast<size_t>(this->out.primitive_matrix) - out_primitive_matrix_i, std::align(types[i].alignment, types[i].length))) );
					out_primitive_matrix_i += std::align(types[i].alignment, types[i].length) * vertices;
				}
			}
		}

		virtual ~primitive() override {
			assert( out.primitive_matrix != nullptr );
				assert( !attributes.empty() );
				::operator delete(out.primitive_matrix, std::align_val_t{std::max(NDCposition.first.alignment,
					__STDCPP_DEFAULT_NEW_ALIGNMENT__)});
		}

		virtual size_t num_arguments() const override { 
			return 1 + attributes.size();
		}

		virtual size_t argument_index(const std::string& name) const override {
			auto attribute_i = attributes.find(name);
			return attribute_i != attributes.end() ? 1 + std::distance(attributes.begin(), attribute_i)
				: name == controlflow::graphics::NDCposition ? 0
				: static_cast<size_t>(-1);
		}

		virtual std::string argument_name(size_t i) const override {
			assert( i < num_results() );
			return i > 0 ? attributes.keys()[i - 1]
				: controlflow::graphics::NDCposition;
		}

		virtual size_t num_results() const override { return num_arguments(); }

		virtual size_t result_index(const std::string& name) const override { return argument_index(name); }

		virtual std::string result_name(size_t i) const override { return argument_name(i); }

		virtual bool set_argument(size_t i, const void* arg) override {
			assert( i < num_arguments() );
			if (i == 0) in.NDCposition = (const unsigned char*)arg;
			else in.attributes[i - 1] = (const unsigned char*)arg;
		}

		virtual const void* get_result(size_t i) const override {
			assert( i < num_results() );
			return i == 0 ? out.primitive_matrix
				: std::next(out.primitive_matrix, attributes.values()[i - 1].second.first/*offset*/);
		}

		virtual void clear() noexcept override {
			index = 0;
		}

		virtual bool setup() override {
			if (vertices == 0) return false;// skip null forloop.
			assert( index == vertices );// no exception when forloop completed.
			assert( setup_vertices != vertices );
			/// move some vertices.
			/// ...
			index = setup_vertices;
			return true;
		}

		virtual void invoke() override {
			std::copy_n(in.NDCposition, NDCposition.first.length,
				std::next(out.primitive_matrix, NDCposition.second * index));
			for (size_t i = 0; i != attributes.size(); ++i) 
				std::copy_n(in.attributes[i], attributes.values()[i].first.length,
					std::next(out.primitive_matrix, attributes.values()[i].second.first + attributes.values()[i].second.second * index));
			++index;
		}

		virtual bool next(size_t) override {
			return index == vertices;
		}
	};

	template<typename vector4>
	struct rfor_each_line : public controlflow::vertex_property {
		using scalar = std::remove_cvref_t<std::declval<vector4>()[0]>;
		struct in_block {
			const vector4* NDCposition;
		} in;
		struct out_block {
			scalar t;
			size_t pixel[2];
		} out;

		size_t xmax, ymax; scalar xmaxf, ymaxf;
		size_t /*x0, y0, */x1, y1; scalar length;
		intptr_t sx, sy, dx, dy, err;

		rfor_each_line(size_t xmax, size_t ymax) : in{nullptr}, out{},
			xmax(xmax), ymax(ymax), xmaxf(static_cast<scalar>(xmax)), ymaxf(static_cast<scalar>(ymax))  {}

		virtual bool setup() override {
			size_t x0 = static_cast<size_t>( std::clamp(round((NDCposition[0][0]*0.5f+1) * xmaxf), scalar(0), xmaxf) );
			size_t y0 = static_cast<size_t>( std::clamp(round((NDCposition[0][1]*0.5f+1) * ymaxf), scalar(0), ymaxf) );
			x1 = static_cast<size_t>( std::clamp(round((NDCposition[1][0]*0.5f+1) * xmaxf), scalar(0), xmaxf) );
			y1 = static_cast<size_t>( std::clamp(round((NDCposition[1][1]*0.5f+1) * ymaxf), scalar(0), ymaxf) );
			if (x0 == x1 && y1 == y1) 
				return false;//skip.

			scalar diffx = NDCposition[1][0] - NDCposition[0][0];
			scalar diffy = NDCposition[1][1] - NDCposition[0][1];
			length = sqrt(diffx*diffx + diffy*diffy);
			dx = x0 < x1 ? intptr_t(x1 - x0) : intptr_t(x0 - x1);
			sx = x0 < x1 ? 1 : -1;
			dy = -(y0 < y1 ? intptr_t(y1 - y0) : intptr_t(y0 - y1));
			sy = y0 < y1 ? 1 : -1;
			err = dx + dy;  /* error value e_xy */
			out.pixel[0] = x0;
			out.pixel[1] = y0;
			return true;
		}

		virtual bool invoke() override {
			scalar diffx = (static_cast<scalar>(out.pixel[0])/xmaxf) * 2 - 1 - NDCposition[0][0];
			scalar diffy = (static_cast<scalar>(out.pixel[1])/ymaxf) * 2 - 1 - NDCposition[0][1];
			out.t = sqrt(diffx*diffx + diffy*diffy)/length;
			if (out.pixel[0] == x1 && out.pixel[1] == y1) {
				return false;//end of loop.
			}

			intptr_t e2 = err * 2;
			if (e2 >= dy) { /* e_xy+e_x > 0 */
				err += dy;
				x0 += sx;
			}
			if (e2 <= dx) { /* e_xy+e_y < 0 */
				err += dx;
				y0 += sy;
			}
			return true;
		}

		virtual void set_argument(size_t i, const void* arg) override {
			assert( i == 0 );
			in.NDCposition = (const vector4*)arg;
		}

		virtual const void* get_result(size_t i) override {
			if (i < sizeof_attributes.size()) {
				auto destination = out.attributes.data();
				for (size_t k = 0; k != i; ++k)
					destination += sizeof_attributes[k];
				return destination;
			} else if (i == sizeof_attributes.size()) {
				return out.pixel;
			}	else {
				abort();
			}
		}
	};

	template<>
	struct rfor_each_line<void> {
		/// ...
	};

	struct write_image_node : public controlflow::vertex_property {
		vec4* image_data;
		size_t image_size[2];
		size_t image_stride_1;
		struct in_block {
			const vec4* color;
			const size_t* pixel;
		} in;

		write_image_node(vec4* image_data, const size_t* image_size) 
			: image_data(image_data), image_size{image_size[0],image_size[1]}, image_stride_1(image_size[0]), in{nullptr,nullptr} {}

		virtual void invoke() override {
			if (in.pixel[0] < image_size[0] && in.pixel[1] < image_size[1]) {
				image_data[in.pixel[0] + in.pixel[1]*image_stride_1] = *in.color;
			}
		}

		virtual void set_argument(size_t i, const void* pointer) override {
			if (i == 0) {
				in.color = (const vec4*)pointer;
			} else if (i == 1) {
				in.pixel = (const size_t*)pointer;
			} else {
				abort();
			}
		}
	};
} }//end of namespace of controlflow::graphics