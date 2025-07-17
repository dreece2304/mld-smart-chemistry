# MLD Smart Chemistry Documentation

## 📚 Documentation Index

### Quick Start
- **[Main README](../README.md)** - Project overview and quick start guide
- **[Installation & Setup](INSTALLATION_SETUP.md)** - Complete installation guide
- **[Usage Guide](USAGE_GUIDE.md)** - How to use the framework

### Remote Work & SSH
- **[SSH & Remote Setup](SSH_REMOTE_SETUP.md)** - Remote execution on lab desktop
- **[SSH Cheat Sheet (HTML)](MLD_SSH_Simulation_Cheatsheet.html)** - Quick reference for SSH commands

### Project Planning & Research
- **[Project Phases](PROJECT_PHASES.md)** - Detailed roadmap and technical specifications
- **[Project Status & Research (HTML)](PROJECT_STATUS_RESEARCH.html)** - Live status dashboard and research priorities

### Technical Documentation
- **[Reorganization Summary](REORGANIZATION_SUMMARY.md)** - Recent project structure changes
- **[API Reference](API_REFERENCE.md)** - *[Coming Soon]* - Code documentation

### Reference Materials
- **[Original Documentation](readmes/)** - Historical documentation for reference

---

## 🔧 How to Update Documentation

### For Project Status & Research (HTML)
The HTML documentation can be updated in two ways:

#### Via Claude
```
Request: "Update the project status to show Phase 3A as 80% complete and add new research priority for X"
```

#### Manual Editing
```html
<!-- Change status indicators -->
<span class="status-indicator status-progress">📝</span>  <!-- In progress -->
<span class="status-indicator status-complete">✓</span>   <!-- Complete -->

<!-- Update progress bars -->
<div class="progress-fill" style="width: 80%"></div>

<!-- Add new research cards -->
<div class="research-card priority-high">
    <h4>🔥 New Research Area</h4>
    <p>Description...</p>
</div>
```

### For Markdown Files
Standard markdown editing for technical documentation:
- Use git for version control
- Keep sections organized with clear headers
- Update timestamps and version numbers

---

## 📖 Documentation Structure

### User Documentation
- **Getting Started**: Installation, setup, first run
- **Usage Guides**: How to use each component
- **Examples**: Common workflows and use cases

### Technical Documentation
- **Architecture**: System design and components
- **API Reference**: Code documentation
- **Development**: Contributing and extending the framework

### Research Documentation
- **Project Phases**: Roadmap and objectives
- **Status Updates**: Current progress and priorities
- **Research Notes**: Findings and future directions

---

## 🎯 Quick Links by Use Case

### New User
1. [Installation & Setup](INSTALLATION_SETUP.md)
2. [Main README](../README.md)
3. [Usage Guide](USAGE_GUIDE.md)

### Remote Work
1. [SSH & Remote Setup](SSH_REMOTE_SETUP.md)
2. [SSH Cheat Sheet (HTML)](MLD_SSH_Simulation_Cheatsheet.html)

### Project Planning
1. [Project Status & Research (HTML)](PROJECT_STATUS_RESEARCH.html)
2. [Project Phases](PROJECT_PHASES.md)

### Development
1. [Reorganization Summary](REORGANIZATION_SUMMARY.md)
2. [API Reference](API_REFERENCE.md) *(Coming Soon)*

---

## 🔄 Documentation Maintenance

### Regular Updates
- **Project Status**: Update progress indicators and research priorities
- **Usage Guide**: Add new features and workflows
- **Installation Guide**: Update for new dependencies

### Version Control
- All documentation is version controlled with git
- Major updates should be committed with descriptive messages
- Use branches for significant documentation restructuring

### HTML vs Markdown
- **HTML**: For dynamic content that needs frequent updates (status, research)
- **Markdown**: For stable technical documentation
- **Both**: Should be kept in sync for overlapping content

---

**Last Updated**: July 17, 2025 | **Version**: 2.0