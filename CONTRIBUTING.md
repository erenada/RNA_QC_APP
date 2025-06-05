# Contributing to RNA Processing APP

Thank you for your interest in contributing to the RNA Processing APP project! This document provides guidelines for contributing to this modular bulk RNA-seq analysis tool suite.

## Getting Started

### Prerequisites
- R version 4.0.0 or higher (4.3.0+ recommended)
- Git
- GitHub account
- Basic knowledge of R, Shiny, and RNA-seq analysis

### Development Setup
1. Fork the repository
2. Clone your fork locally
3. Set up the development environment:
   ```r
   # Navigate to the project directory
   setwd("path/to/APP")
   
   # Install dependencies (automatic on first app launch)
   shiny::runApp("app.R")
   ```

## How to Contribute

### Reporting Issues
1. Check existing issues first
2. Use the appropriate issue template
3. Include:
   - System information (R version, OS, browser)
   - Reproducible example
   - Expected vs. actual behavior
   - Screenshots if applicable

### Suggesting Features
1. Open an issue with the "Feature Request" label
2. Clearly describe the feature and its benefits
3. Consider the modular architecture when proposing new features

### Code Contributions

#### Branch Naming
- `feature/module-name-description` for new features
- `fix/issue-description` for bug fixes
- `docs/documentation-updates` for documentation changes

#### Code Standards
- Follow R style guidelines (use `styler` package)
- Include comprehensive comments
- Write unit tests for new functions
- Update documentation

#### Pull Request Process
1. Create a feature branch from `main`
2. Make your changes
3. Add tests for new functionality
4. Update documentation
5. Ensure all tests pass
6. Submit pull request with:
   - Clear description of changes
   - Reference to related issues
   - Screenshots of UI changes (if applicable)

### Module Development Guidelines

#### Architecture Principles
- Each module should be self-contained
- Use consistent naming conventions
- Follow the established UI/Server pattern
- Include comprehensive error handling

#### File Structure
```
modules/
├── module_name/
│   ├── R/
│   │   ├── ui_*.R
│   │   ├── server_*.R
│   │   └── utils_*.R
│   ├── tests/
│   │   └── testthat/
│   └── www/
│       └── styles.css
```

#### Testing Requirements
- Write unit tests for all utility functions
- Include integration tests for module functionality
- Test with example data
- Verify cross-browser compatibility

## Documentation Standards

### Code Documentation
- Use roxygen2 for function documentation
- Include examples in function documentation
- Document all parameters and return values

### User Documentation
- Update user manual for new features
- Include screenshots for UI changes
- Provide example workflows
- Update troubleshooting guide if needed

### Developer Documentation
- Document architectural decisions
- Update API documentation
- Include development roadmap updates

## Quality Assurance

### Before Submitting
- [ ] Code follows style guidelines
- [ ] All tests pass
- [ ] Documentation is updated
- [ ] No new linter warnings
- [ ] UI is responsive and accessible
- [ ] Performance impact is considered

### Review Process
1. Automated checks run on PR
2. Code review by maintainers
3. Testing in different environments
4. Documentation review
5. Final approval and merge

## Community Guidelines

### Communication
- Be respectful and constructive
- Use clear, professional language
- Focus on technical merit
- Help others learn and improve

### Recognition
Contributors will be acknowledged in:
- CONTRIBUTORS.md file
- Release notes
- Documentation credits

## Getting Help

### Resources
- [Developer Documentation](docs/developer/README.md)
- [User Manual](docs/user_manual/README.md)
- [Project Status](docs/project_status.md)

### Contact
- Create an issue for bugs or feature requests
- Use discussions for general questions
- Email: [maintainer email if applicable]

## License

By contributing to this project, you agree that your contributions will be licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.

---

**Author:** Eren Ada, PhD  
**Last Updated:** 6/3/2025 