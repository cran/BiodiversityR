.onAttach <- function(...) {
      packageStartupMessage("BiodiversityR ", utils::packageDescription("BiodiversityR", field="Version"),
      ": Use command BiodiversityRGUI() to launch the Graphical User Interface; \nto see changes use BiodiversityRGUI(changeLog=TRUE, backward.compatibility.messages=TRUE)\n")
}

