`BiodiversityR.changeLog` <-
function()
{
    change.file <- file.path(system.file(package = "BiodiversityR"), "ChangeLog")
    file.show(change.file)
}

