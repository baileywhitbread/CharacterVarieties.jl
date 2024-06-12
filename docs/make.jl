using CharacterVarieties
using Documenter

DocMeta.setdocmeta!(CharacterVarieties, :DocTestSetup, :(using CharacterVarieties); recursive=true)

makedocs(;
    modules=[CharacterVarieties],
    authors="Bailey Whitbread",
    sitename="CharacterVarieties.jl",
    format=Documenter.HTML(;
        canonical="https://baileywhitbread.github.io/CharacterVarieties.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/baileywhitbread/CharacterVarieties.jl",
    devbranch="master",
)
