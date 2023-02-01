#!/bin/sh
black . -v
flake8 . --count --ignore="E501" --show-source --statistics --exclude=./.*
