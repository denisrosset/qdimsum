name := "qdimsum"

// orgnization name (e.g., the package name of the project)
organization := "com.faacets"

version := "0.1-SNAPSHOT"

// Enables publishing to maven repo
publishMavenStyle := true

// Do not append Scala versions to the generated artifacts
crossPaths := false

// This forbids including Scala related libraries into the dependency
autoScalaLibrary := false

javacOptions ++= Seq("-source", "1.6", "-target", "1.6")
