
using GenieFramework
using HTTP
include("app.jl")
up(async=true)
HTTP.get("http://localhost:8000/");
