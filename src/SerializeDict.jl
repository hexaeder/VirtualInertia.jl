using Serialization

export SerializeDict

struct SerializeDict
    keys::Set{UInt}
    dir::String
end

_keypath(sd::SerializeDict, k) = joinpath(sd.dir, string(hash(k)))

function SerializeDict()
    dir = tempname()
    mkdir(dir)
    # @info "Start SerializeDict at $dir"
    SerializeDict(Set{Any}(), dir)
end

function SerializeDict(dir)
    isdir(dir) || mkdir(dir)
    SerializeDict(Set{Any}(), dir)
end

Base.haskey(sd::SerializeDict, k) = hash(k) ∈ sd.keys

function Base.getindex(sd::SerializeDict, k)
    if haskey(sd, k)
        path = _keypath(sd, k)
        # @info "Load cached" k
        return deserialize(path)
    else
        throw(KeyError(k))
    end
end

function Base.setindex!(sd::SerializeDict, v, k)
    # @info "Pushed key" k
    hk = hash(k)
    path = _keypath(sd, k)
    serialize(path, v)
    push!(sd.keys, hk)
    v
end

function Base.delete!(sd::SerializeDict, k)
    path = _keypath(sd, k)
    isfile(path) && rm(path)
    delete!(sd.keys, hash(k))
    sd
end

function Base.empty!(sd::SerializeDict)
    # @info "Dict emptied"
    for hk in sd.keys
        path = joinpath(sd.dir, string(hk))
        rm(path)
    end
    empty!(sd.keys)
    sd
end

Base.isempty(sd::SerializeDict) = isempty(sd.keys)
Base.length(sd::SerializeDict) = length(sd.keys)

function Base.get(sd::SerializeDict, key, default)
    haskey(sd, key) ? sd[key] : default
end

function Base.get!(sd::SerializeDict, key, default)
    if haskey(sd, key)
        return sd[key]
    else
        sd[key] = default
        return default
    end
end

function Base.get!(f::Union{Function, Type}, sd::SerializeDict, key)
    if haskey(sd, key)
        return sd[key]
    else
        val = f()
        sd[key] = val
        return val
    end
end

# function Base.get(sd::SerializeDict, k, default)
# end

# function Base.getkey
# Base.getindex, Base.get, Base.haskey, Base.getkey, Base.isempty, Base.length
# Base.delete!, Base.empty!, Base.setindex!
