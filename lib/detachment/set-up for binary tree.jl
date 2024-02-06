struct HeapElement{T}
    value::T
    coordinates::Tuple{Int, Int}
end

function Base.isless(x::HeapElement, y::HeapElement)
    return isless(x.value, y.value)
end

export HeapElement

using DataStructures
using .Struct_Module
# Create a binary heap of Struct_Module.HeapElement

a = Struct_Module.HeapElement(5, (3, 9))
b = Struct_Module.HeapElement(10, (8, 9))
a < b

