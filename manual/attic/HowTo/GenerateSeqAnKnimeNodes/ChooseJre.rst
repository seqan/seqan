In order to change the JRE to be used by KNIME go to Eclipse
MenuTrace(Preferences) and select the MenuTrace(Java) menu.

`Image(Preferences.png, width=400,
align=center) <Image(Preferences.png, width=400, align=center)>`__

Afterwards you can add the right JRE. Under MacOs you choose the entry
`MenuTrace(MacOS X VM) <MenuTrace(MacOS X VM)>`__

`Image(MacOSJre.png, width=400,
align=center) <Image(MacOSJre.png, width=400, align=center)>`__

then press next and select the right path, which should be
*/System/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Home*

ash shown here:

`Image(JreHome.png, width=600,
align=center) <Image(JreHome.png, width=600, align=center)>`__

Press MenuTrace(Finish) and the right JRE will be added.

Afterwards you have to set the compiler options. In order to do so go to
MenuTrace(Project) and select MenuTrace(Properties)

`Image(ProjectProperties.png, width=400,
align=center) <Image(ProjectProperties.png, width=400, align=center)>`__

No select `MenuTrace(Java Compiler) <MenuTrace(Java Compiler)>`__ and
select the correct JRE at the `MenuTrace(Compiler compliance
level:) <MenuTrace(Compiler compliance level:)>`__

`Image(JavaCompiler.png, width=400,
align=center) <Image(JavaCompiler.png, width=400, align=center)>`__

If you run the project now KNIME should start without problems.

.. raw:: mediawiki

   {{TracNotice|{{PAGENAME}}}}
