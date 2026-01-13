 find . -type f \( -name "*.py" -o -name "*.sh" -o -name "./example/*.yaml" \) | npx repomix@latest --stdin
