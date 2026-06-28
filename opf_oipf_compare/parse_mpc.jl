# -----------------------------------------------------------------------------
# Lightweight MATPOWER .m case parser for Julia.
# -----------------------------------------------------------------------------

function _strip_matlab_comments(line::String)
    idx = findfirst('%', line)
    return isnothing(idx) ? line : line[1:prevind(line, idx)]
end

function _strip_matlab_comments(text::AbstractString)
    lines = split(String(text), '\n')
    cleaned = String[]

    for line in lines
        idx = findfirst('%', line)
        if idx === nothing
            push!(cleaned, line)
        else
            push!(cleaned, line[1:prevind(line, idx)])
        end
    end

    return join(cleaned, "\n")
end

function _parse_matrix(content::String, field::String)
    pattern = Regex("mpc\\.$field\\s*=\\s*\\[(.*?)\\];", "s")
    m = match(pattern, content)
    if isnothing(m)
        return zeros(Float64, 0, 0)
    end

    block = m.captures[1]
    rows = Vector{Vector{Float64}}()

    for raw_line in split(block, '\n')
        line = _strip_matlab_comments(raw_line)
        line = strip(replace(line, ";" => " ", "," => " "))

        if isempty(line)
            continue
        end

        vals = split(line)
        if isempty(vals)
            continue
        end

        push!(rows, parse.(Float64, vals))
    end

    if isempty(rows)
        return zeros(Float64, 0, 0)
    end

    ncol = length(rows[1])
    mat = zeros(Float64, length(rows), ncol)

    for i in eachindex(rows)
        if length(rows[i]) != ncol
            error("Inconsistent number of columns in mpc.$field at row $i.")
        end
        mat[i, :] .= rows[i]
    end

    return mat
end

function parse_mpc(path::AbstractString)
    content = read(path, String)

    mpc = Dict{String, Any}()

    version = _parse_scalar(content, "version")
    baseMVA = _parse_scalar(content, "baseMVA")

    if !isnothing(version)
        mpc["version"] = version
    end

    if isnothing(baseMVA)
        error("Cannot find mpc.baseMVA in file: $path")
    end

    mpc["baseMVA"] = baseMVA
    mpc["bus"]     = _parse_matrix(content, "bus")
    mpc["gen"]     = _parse_matrix(content, "gen")
    mpc["branch"]  = _parse_matrix(content, "branch")
    mpc["gencost"] = _parse_matrix(content, "gencost")

    return mpc
end

function _strip_matlab_comments(text::AbstractString)
    lines = split(String(text), '\n')
    cleaned = String[]

    for line in lines
        idx = findfirst('%', line)
        push!(cleaned, idx === nothing ? line : line[1:prevind(line, idx)])
    end

    return join(cleaned, "\n")
end

function _parse_scalar(content::AbstractString, field::String)
    pattern = Regex("mpc\\.$field\\s*=\\s*([^;]+);")
    m = match(pattern, String(content))

    if isnothing(m)
        return nothing
    end

    value = strip(m.captures[1])
    value = replace(value, "'" => "", "\"" => "")

    if field == "version"
        return value
    else
        return parse(Float64, value)
    end
end

function _parse_matrix(content::AbstractString, field::String)
    content_clean = _strip_matlab_comments(content)

    pattern = Regex("mpc\\.$field\\s*=\\s*\\[(.*?)\\];", "s")
    m = match(pattern, content_clean)

    if isnothing(m)
        return zeros(Float64, 0, 0)
    end

    block = m.captures[1]
    rows = Vector{Vector{Float64}}()

    for raw_line in split(block, '\n')
        line = strip(replace(raw_line, ";" => " ", "," => " "))

        if isempty(line)
            continue
        end

        vals = split(line)

        if isempty(vals)
            continue
        end

        push!(rows, parse.(Float64, vals))
    end

    if isempty(rows)
        return zeros(Float64, 0, 0)
    end

    ncol = length(rows[1])
    mat = zeros(Float64, length(rows), ncol)

    for i in eachindex(rows)
        if length(rows[i]) != ncol
            error("Inconsistent number of columns in mpc.$field at row $i.")
        end
        mat[i, :] .= rows[i]
    end

    return mat
end

function parse_mpc(path::AbstractString)
    content = read(path, String)

    mpc = Dict{String, Any}()

    version = _parse_scalar(content, "version")
    baseMVA = _parse_scalar(content, "baseMVA")

    if !isnothing(version)
        mpc["version"] = version
    end

    if isnothing(baseMVA)
        error("Cannot find mpc.baseMVA in file: $path")
    end

    mpc["baseMVA"] = baseMVA
    mpc["bus"]     = _parse_matrix(content, "bus")
    mpc["gen"]     = _parse_matrix(content, "gen")
    mpc["branch"]  = _parse_matrix(content, "branch")
    mpc["gencost"] = _parse_matrix(content, "gencost")

    return mpc
end