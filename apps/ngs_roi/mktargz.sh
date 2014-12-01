#!/bin/bash

pushd tool_shed
tar --exclude=.svn --exclude="*.pyc" --exclude="#*" --exclude=".*" --exclude="*~" -z -c -v -f ../ngs_roi.tar.gz *
popd
