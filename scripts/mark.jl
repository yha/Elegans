#using Revise
using Elegans, Atom
using Interact: textbox
#Atom.displayinplotpane(textbox("")) && sleep(1) # workaround Interact.jl issue #321
mark_stages_gui()

##

using Elegans
mark_stages_window()
