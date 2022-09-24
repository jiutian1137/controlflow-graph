
template<typename _Ty, size_t _Size>
struct vecN { 
	_Ty data[_Size];
	size_t size() const { return _Size; }
	_Ty operator[](size_t i) const { return data[i]; }
	_Ty& operator[](size_t i) { return data[i]; }
};

template<typename _Ty, size_t _Size>
inline _Ty dot(const vecN<_Ty,_Size>& x, const vecN<_Ty,_Size>& y) { 
	_Ty _Dot = x[0] * y[0];
	for (size_t i = 1; i != _Size; ++i)
		_Dot += x[i] * y[i];
	return _Dot;
}

template<typename _Ty, size_t _Size>
inline _Ty length(const vecN<_Ty,_Size>& x) { 
	_Ty _Dot = x[0] * x[0];
	for (size_t i = 1; i != _Size; ++i)
		_Dot += x[i] * x[i];
	return sqrt(_Dot);
}

template<typename _Ty, size_t _Size>
inline vecN<_Ty,_Size> operator+(const vecN<_Ty,_Size>& x, const vecN<_Ty,_Size>& y) {
	vecN<_Ty,_Size> z;
	for (size_t i = 0; i != _Size; ++i)
		z[i] = x[i] + y[i];
	return z;
}

template<typename _Ty, size_t _Size>
inline vecN<_Ty,_Size> operator-(const vecN<_Ty,_Size>& x, const vecN<_Ty,_Size>& y) {
	vecN<_Ty,_Size> z;
	for (size_t i = 0; i != _Size; ++i)
		z[i] = x[i] - y[i];
	return z;
}

template<typename _Ty, size_t _Size>
inline vecN<_Ty,_Size> operator*(const vecN<_Ty,_Size>& x, const vecN<_Ty,_Size>& y) {
	vecN<_Ty,_Size> z;
	for (size_t i = 0; i != _Size; ++i)
		z[i] = x[i] * y[i];
	return z;
}

template<typename _Ty, size_t _Size>
inline vecN<_Ty,_Size> operator*(const vecN<_Ty,_Size>& x, const _Ty& y) {
	vecN<_Ty,_Size> z;
	for (size_t i = 0; i != _Size; ++i)
		z[i] = x[i] * y;
	return z;
}

template<typename _Ty, size_t _Size>
inline vecN<_Ty,_Size> operator/(const vecN<_Ty,_Size>& x, const vecN<_Ty,_Size>& y) {
	vecN<_Ty,_Size> z;
	for (size_t i = 0; i != _Size; ++i)
		z[i] = x[i] / y[i];
	return z;
}

using vec2 = vecN<float, 2>;
using vec3 = vecN<float, 3>;
using vec4 = vecN<float, 4>;

template<typename _Ty, typename _Ty2>
_Ty lerp(const _Ty& x, const _Ty& y, const _Ty2 t) {
	return (y - x) * t + x;
}

#if 0
struct constants_node : public cfg_vertex_property {
	std::string names;
	std::vector<size_t> sizeof_properties;
	std::vector<size_t> offsetof_properties;
	struct out_block {
		std::vector<char> buffer;
	} out;

	virtual const void* get_result(size_t i) override {
		return &out.buffer[offsetof_properties[i]];
	}
};

template<typename _Fn, typename _Ty1, typename _Ty2 = _Ty1>
struct binaryfn_node : public cfg_vertex_property {
	_Fn binaryfn;
	struct in_block {
		const _Ty1* left;
		const _Ty2* right;
	} in;
	struct out_block {
		decltype(std::declval<_Fn>()(std::declval<_Ty1>(), std::declval<_Ty2>())) result;
	} out;

	explicit binaryfn_node(const _Fn& binaryfn) : binaryfn(binaryfn), in{ nullptr,nullptr }, out{} {}

	virtual void invoke() override {
		out.result = binaryfn(*in.left, *in.right);
	}

	virtual void set_argument(size_t i, const void* pointer) override {
		if (i == 0) {
			in.left = (const _Ty1*)pointer;
		} else if (i == 1) {
			in.right = (const _Ty2*)pointer;
		} else {
			abort();
		}
	}

	virtual const void* get_result(size_t i) override {
		assert( i == 0 );
		return &out.result;
	}
};

template<typename _Ty>
struct constant_node : public cfg_vertex_property {
	struct out_block {
		_Ty value;
	} out;

	constexpr constant_node(const _Ty& value) : out{ value } {}

	constexpr constant_node(_Ty&& value) : out{ std::move(value) } {}

	virtual const void* get_result(size_t i) override {
		assert( i == 0 );
		return &out.value;
	}
};

struct opengl_glsl_text {
	std::string global;
	std::string local;
};
#endif


#include <iostream>

#include "rasterization.hpp"
#include "geometry/graph.hpp"
#include "controlflow/basic.hpp"

#include <functional>

int main() {
	/*for (size_t i = 2; i != 100; ++i) {
		bool prime = true;
		for (size_t j = 2; j < i; ++j) {
			prime &= (i % j != 0);
		}

		if (prime) {
			std::cout << i << std::endl;
		}
	}

	return 0;*/

	geometry::adjacency_list<geometry::bidirected, false, true, true,
		std::string, std::shared_ptr<controlflow::vertex_property>, controlflow::edge_property> g;

	g.add_vertex("Print Int", std::make_shared<controlflow::print<int>>(std::_Ptr_cout, ' '));

	g.add_vertex("Variables", std::make_shared<controlflow::buffer>(std::vector<std::string>{
		"count", "zero", "newline", "is prime", "true" },
		 size_t(100), size_t(0), char('\n'), bool(false), bool(true) )
	);
	
	g.add_vertex("outer for", std::make_shared<controlflow::for_>(2,1000));
	g.add_edge(g.find("Variables"), g.find("outer for"), { {}, true });

	g.add_vertex("clear", std::make_shared<controlflow::assign<bool>>());
	g.add_edge(g.find("Variables"), g.find("clear"), {{ 
		{ g[g.find("Variables")].prop->result_index("is prime"), g[g.find("clear")].prop->argument_index("left") },
		{ g[g.find("Variables")].prop->result_index("true"), g[g.find("clear")].prop->argument_index("right") }
		}, false});
	g.add_edge(g.find("outer for"), g.find("clear"), { {}, true });

	g.add_vertex("inner for", std::make_shared<controlflow::for_>(2));
	g.add_edge(g.find("outer for"), g.find("inner for"), {{ 
		{ g[g.find("outer for")].prop->result_index("index"), g[g.find("inner for")].prop->argument_index("end") } }, true });
	
#if 0
	g.add_edge(g.find("inner for"), g.find("Print Int"), {{ 
		{ g[g.find("inner for")].prop->result_index("index"), 0 } },true});

	g.add_vertex("Print Char", std::make_shared<controlflow::print<char>>(std::_Ptr_cout, ' '));
	g.add_edge(g.find("Variables"), g.find("Print Char"), { { { g[g.find("Variables")].prop->result_index("newline"), 0 } }, false });
	g.add_edge(g.find("outer for"), g.find("Print Char"), { {},true });
#endif

	g.add_vertex("modulus", std::make_shared<controlflow::modulus<size_t>>());
	g.add_edge(g.find("outer for"), g.find("modulus"), {{ 
		{ g[g.find("outer for")].prop->result_index("index"), g[g.find("modulus")].prop->argument_index("left") } }, false});
	g.add_edge(g.find("inner for"), g.find("modulus"), {{ 
		{ g[g.find("inner for")].prop->result_index("index"), g[g.find("modulus")].prop->argument_index("right") } }, true});
	
	g.add_vertex("is prime", std::make_shared<controlflow::not_equal_to<size_t>>());
	g.add_edge(g.find("modulus"), g.find("is prime"), {{ 
		{ g[g.find("modulus")].prop->result_index("result"), g[g.find("is prime")].prop->argument_index("left") } }, true});
	g.add_edge(g.find("Variables"), g.find("is prime"), {{ 
		{ g[g.find("Variables")].prop->result_index("zero"), g[g.find("is prime")].prop->argument_index("right") } }, false});

	g.add_vertex("accumulate", std::make_shared<controlflow::assign_bit_and<bool>>());
	g.add_edge(g.find("Variables"), g.find("accumulate"), {{ 
		{ g[g.find("Variables")].prop->result_index("is prime"), g[g.find("accumulate")].prop->argument_index("left") } }, false});
	g.add_edge(g.find("is prime"), g.find("accumulate"), {{ 
		{ g[g.find("is prime")].prop->result_index("result"), g[g.find("accumulate")].prop->argument_index("right") } }, true});

	g.add_vertex("check prime", std::make_shared<controlflow::if_>());
	g.add_edge(g.find("Variables"), g.find("check prime"), { { { g[g.find("Variables")].prop->result_index("is prime"), 0 } }, false });
	g.add_edge(g.find("outer for"), g.find("check prime"), { {}, true });

	g.add_vertex("will do check prime", std::make_shared<controlflow::case_>());
	g.add_edge(g.find("check prime"), g.find("will do check prime"), { {
		{ 0, 0 } },true });

	g.add_edge(g.find("will do check prime"), g.find("Print Int"), { {}, true });
	g.add_edge(g.find("outer for"), g.find("Print Int"), { {
		{ g[g.find("outer for")].prop->result_index("index"), 0 } },false });

#if 0
	vec4 image[64][64] = { 0 };
	size_t image_size[2] = { 64,64 };
	size_t image_edge[2] = { 63,63 };

	std::vector<vec4> vertices = {
		/// Lines[0]
		vec4{  0,  0,0,1},
		vec4{0.5,0.5,0,1},
		/// Lines[1]
		vec4{0.5,0.5,0,1},
		vec4{1.0,0.5,0,1},
	};

	g.add_vertex("Vertex Node", 
		std::make_shared<for_each_vertex_n>(std::initializer_list{ attribute_view{(const char*)vertices.data(),sizeof(vertices[0]),sizeof(vertices[0])} }, vertices.size())
		);
	g.add_vertex("Primitive Node",
		std::make_shared<primitive>(2, std::initializer_list{ sizeof(vec4) })
		);
	g.add_vertex("Raster Node",
		std::make_shared<raster_primitive>(2, std::initializer_list{ sizeof(vec4) }, image_edge)
		);
	g.add_vertex("Output",
		std::make_shared<write_image_node>((vec4*)image, image_size)
		);

	g.add_edge(g.find("Vertex Node"), g.find("Primitive Node"), 
		{{ 0/*attributes[0]*/, 0/*primitive.attributes[0]*/ }});
	g.add_edge(g.find("Primitive Node"), g.find("Raster Node"), 
		{{ 0/*primitive.attributes[0]*/, 0/*attributes[0]*/ }});

	g.add_edge(g.find("Raster Node"), g.find("Output"), 
		{{ 1/*pixel*/, 1/*pixel*/ }});

	g.add_vertex("Exposure", std::make_shared<constant_node<float>>(10.0f));
	g.add_vertex("Exposured Color", std::make_shared<binaryfn_node<vec4(*)(vec4,float), vec4, float>>([](vec4 a, float s){ return a*s; }));

	g.add_edge(g.find("Raster Node"), g.find("Exposured Color"), 
		{{ 0/*attributes[0]*/, 0/*left*/}});
	g.add_edge(g.find("Exposure"), g.find("Exposured Color"),
		{{ 0/*value*/, 1/*right*/ }});
	g.add_edge(g.find("Exposured Color"), g.find("Output"), 
		{{ 0/*result*/, 0/*color*/ }});

	///@note link don't cross seperate stages (exmaple: "Vertex Node" -> "Exposure" is ERROR). 
	g.add_edge(g.find("Raster Node"), g.find("Exposure"));
#endif

	/// set_argument.
	std::queue<decltype(g)::vertex_descriptor> Q;
	for (auto vfirst = g.c.vertices.begin(); vfirst != g.c.vertices.end(); ++vfirst) 
		if ((*vfirst).second.in_edges.empty()) 
			Q.push(g._Get_descriptor(g.c.vertices, vfirst));
	//assert(Q.size() == 1);
	while (!Q.empty()) {
		auto source = Q.front();
		for (auto eiter = g[source].out_edges.begin(); eiter != g[source].out_edges.end(); ++eiter) {
			auto e = g[g.find(source, (*eiter).first)];
			for (const auto& link : e.prop.first) {
				e.target.prop->set_argument(link.second, 
					e.source.prop->get_result(link.first));
			}
			Q.push((*eiter).first);
		}

		Q.pop();
	}

	/// run.
	//std::map<decltype(g)::vertex_descriptor,  std::set<decltype(g)::vertex_descriptor>> lock;
	std::stack<decltype(g)::vertex_descriptor> S;
	for (auto v = g.c.vertices.begin(); v != g.c.vertices.end(); ++v) 
		if ((*v).second.in_edges.empty()) 
			S.push(g._Get_descriptor(g.c.vertices, v));
	//assert(S.size() == 1);
	while (!S.empty()) {
		auto source = S.top();
		
		controlflow::invoke_result_t result = g[source].prop->invoke();
		
		if (result & controlflow::invoke_throw) {
			/// ...
		}
		
		if (result & controlflow::invoke_pop_current) {
			S.pop();
		}
		
		if (result & controlflow::invoke_pop_last) {
			assert(!S.empty());
			S.pop();
		}

		if (result & controlflow::invoke_next) {
			if (result & controlflow::invoke_nextpos) {
				size_t source_slot = (result>>controlflow::invoke_nextpos_shift);
				for (auto e = g[source].out_edges.rbegin(); e != g[source].out_edges.rend(); ++e) {
					assert( (*e).second->second );
					if ( std::any_of((*e).second->first.begin(), (*e).second->first.end(), [source_slot](const std::pair<size_t,size_t>& link){ return link.first == source_slot; }) ) {
						auto target = (*e).first;
						if ( g[target].prop->setup() ) {
							S.push(target);
						}
					}
				}
			} else {
				for (auto e = g[source].out_edges.rbegin(); e != g[source].out_edges.rend(); ++e) {
					auto target = (*e).first;
					if ((*e).second->second) {
						if ( g[target].prop->setup() ) {
							S.push(target);
						}
					}
				}
			}
		}

#if 0
		/// example forloop, 
		/// when the node not completed entire loop (i != N),
		/// the node should remain at original position (at front), that it can be invoked until completed entire loop.
		/// so we cannot use breadth-first-search, because the queue will be blocked (front cannot be pop).
		if (source.prop->pop()) 
			S.pop();
#endif

#if 0
		for (auto e = source.out_edges.rbegin(); e != source.out_edges.rend(); ++e) 
			if ( ((*e).second->empty() && source.prop->next(-1))
				|| std::any_of((*e).second->begin(), (*e).second->end(), [&source](const std::pair<size_t, size_t>& link) { return source.prop->next(link.first); }))
				if (++lock[(*e).first] == g[(*e).first].in_edges.size()) {
					lock[(*e).first] = 0;//reset
					if ( g[(*e).first].prop->setup() ) {
						S.push((*e).first);
					}
				}
#endif
	}

#if 0
	geometry::adjacency_list<geometry::bidirected, true, false, false, int, int> g10;
	geometry::adjacency_list<geometry::bidirected, true, true, false, int, int> g11;
	geometry::adjacency_list<geometry::bidirected, true, true, true, int, int> g12;
	geometry::adjacency_list<geometry::bidirected, false, true, false, int, int> g20;
	geometry::adjacency_list<geometry::bidirected, false, true, true, int, int> g21;
	geometry::adjacency_list<geometry::bidirected, false, false, true, int, int> g30;
	geometry::adjacency_list<geometry::bidirected, true, false, true, int, int> g31;

	geometry::adjacency_list<geometry::bidirected, true, false, false, std::string, int> lg10;
	geometry::adjacency_list<geometry::bidirected, true, true, false, std::string, int> lg11;
	geometry::adjacency_list<geometry::bidirected, true, true, true, std::string, int> lg12;
	geometry::adjacency_list<geometry::bidirected, false, true, false, std::string, int> lg20;
	geometry::adjacency_list<geometry::bidirected, false, true, true, std::string, int> lg21;
	geometry::adjacency_list<geometry::bidirected, false, false, true, std::string, int> lg30;
	geometry::adjacency_list<geometry::bidirected, true, false, true, std::string, int> lg31;

	std::cout << "Test unlabeled_graph::add_vertex(prop)    {\n";
	{
		std::cout << g10.contains(10) << ',' << g11.contains(10) << ',' << g12.contains(10) << std::endl;
		std::cout << g20.contains(10) << ',' << g21.contains(10) << std::endl;
		std::cout << g30.contains(10) << ',' << g31.contains(10) << std::endl;

		for (int i = 0; i != 30; ++i) {
			g10.add_vertex(i); g11.add_vertex(i); g12.add_vertex(i);
			g20.add_vertex(i); g21.add_vertex(i);
			g30.add_vertex(i); g31.add_vertex(i);
		}

		std::cout << g10.contains(10) << ',' << g11.contains(10) << ',' << g12.contains(10) << std::endl;
		std::cout << g20.contains(10) << ',' << g21.contains(10) << std::endl;
		std::cout << g30.contains(10) << ',' << g31.contains(10) << std::endl;
	}
	std::cout << "}\n";
	
	std::cout << "Test labeled_graph::add_vertex(id, prop)    {\n";
	{
		std::cout << lg10.contains("10") << ',' << lg11.contains("10") << ',' << lg12.contains("10") << std::endl;
		std::cout << lg20.contains("10") << ',' << lg21.contains("10") << std::endl;
		std::cout << lg30.contains("10") << ',' << lg31.contains("10") << std::endl;

		for (int i = 0; i != 30; ++i) {
			lg10.add_vertex(std::to_string(i), rand()); lg11.add_vertex(std::to_string(i), rand()); lg12.add_vertex(std::to_string(i), rand());
			lg20.add_vertex(std::to_string(i), rand()); lg21.add_vertex(std::to_string(i), rand());
			lg30.add_vertex(std::to_string(i), rand()); lg31.add_vertex(std::to_string(i), rand());
		}

		std::cout << lg10.contains("10") << ',' << lg11.contains("10") << ',' << lg12.contains("10") << std::endl;
		std::cout << lg20.contains("10") << ',' << lg21.contains("10") << std::endl;
		std::cout << lg30.contains("10") << ',' << lg31.contains("10") << std::endl;

		auto sss = lg10.add_edge(lg10.find("10"), lg10.find("20"));
		lg10.remove_vertex(lg10.find("20"));
		//lg10.remove_edge(lg10.find("10"), lg10.find("20"));
		std::cout << "done.\n";
	}
	std::cout << "}\n";
#endif

	return 0;
}