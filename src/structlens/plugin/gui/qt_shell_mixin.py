"""Focused shell presentation/routing mixin for the Qt panel."""

# ruff: noqa: F403,F405

from __future__ import annotations

from structlens import __version__

from .model import PAGE_LABELS
from .qt_context import *  # noqa: F401,F403


class ShellMixin(QtMixinContext):
    """Cohesive GUI-only shell behavior composed into PanelController."""

    def _build_shell(self) -> None:
        self.widget.setStyleSheet(_stylesheet())
        root = self.w.QVBoxLayout(self.widget)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        header = self.w.QFrame(self.widget)
        header.setObjectName("header")
        header_layout = self.w.QHBoxLayout(header)
        header_layout.setContentsMargins(16, 12, 16, 12)
        header_layout.setSpacing(12)
        brand = self.w.QVBoxLayout()
        brand.setSpacing(2)
        brand.addWidget(_label(self.w, "StructLens", "windowTitle"))
        self.workflow_context = _label(self.w, "Choose reference and target structures to begin.", "workflowContext")
        self.workflow_context.setWordWrap(True)
        size_policy = getattr(self.w.QSizePolicy, "Policy", self.w.QSizePolicy)
        self.workflow_context.setSizePolicy(size_policy.Ignored, size_policy.Preferred)
        brand.addWidget(self.workflow_context)
        header_layout.addLayout(brand, 1)
        self.header_status = _label(self.w, "Ready", "statusPill")
        header_layout.addWidget(self.header_status)
        self.cancel_button = _button(self.w, "Cancel", "secondaryButton")
        self.cancel_button.setVisible(False)
        self.cancel_button.clicked.connect(self._cancel_analysis)
        header_layout.addWidget(self.cancel_button)
        self.compare_button = _button(self.w, "Compare structures", "primaryButton")
        header_layout.addWidget(self.compare_button)
        root.addWidget(header)

        body = self.w.QWidget(self.widget)
        body_layout = self.w.QHBoxLayout(body)
        body_layout.setContentsMargins(0, 0, 0, 0)
        body_layout.setSpacing(0)
        self.sidebar = self.w.QFrame(body)
        self.sidebar.setObjectName("sidebar")
        self.sidebar.setFixedWidth(178)
        side_layout = self.w.QVBoxLayout(self.sidebar)
        side_layout.setContentsMargins(10, 14, 10, 14)
        side_layout.setSpacing(8)
        side_layout.addWidget(_label(self.w, "WORKFLOW", "sectionKicker"))
        self.nav = self.w.QListWidget(self.sidebar)
        self.nav.setObjectName("workflowNav")
        self.nav.setSpacing(2)
        for section in SCIENTIFIC_SECTIONS:
            item = self.w.QListWidgetItem("2. Configure" if section == "Structures" else PAGE_LABELS.get(section, section))
            item.setToolTip(_page_subtitle(section, standalone=self.command is None))
            self.nav.addItem(item)
        self.nav.setCurrentRow(0)
        side_layout.addWidget(self.nav, 1)
        sidebar_note = _label(self.w, "Load → Configure → Compare\nReview → Export", "sidebarNote")
        sidebar_note.setWordWrap(True)
        side_layout.addWidget(sidebar_note)
        body_layout.addWidget(self.sidebar)

        self.pages = self.w.QStackedWidget(body)
        self.pages.setObjectName("pageStack")
        body_layout.addWidget(self.pages, 1)
        root.addWidget(body, 1)

        footer = self.w.QFrame(self.widget)
        footer.setObjectName("footer")
        footer_layout = self.w.QHBoxLayout(footer)
        footer_layout.setContentsMargins(20, 8, 20, 8)
        self.footer_status = _label(self.w, "Choose two structures to begin.", "footerStatus")
        self.footer_status.setWordWrap(True)
        self.footer_status.setSizePolicy(size_policy.Ignored, size_policy.Preferred)
        footer_layout.addWidget(self.footer_status, 1)
        self.progress = self.w.QProgressBar(footer)
        self.progress.setObjectName("analysisProgress")
        self.progress.setRange(0, 0)
        self.progress.setMaximumWidth(120)
        self.progress.setVisible(False)
        footer_layout.addWidget(self.progress)
        self.footer_meta = _label(self.w, f"v{__version__}", "footerMeta")
        footer_layout.addWidget(self.footer_meta)
        root.addWidget(footer)

        self.compare_button.clicked.connect(self._start_analysis)

    def _add_page(self, title: str, purpose: str) -> tuple[Any, Any]:
        page = self.w.QWidget(self.pages)
        page.setObjectName(f"page{title}")
        outer = self.w.QVBoxLayout(page)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        scroll = self.w.QScrollArea(page)
        scroll.setObjectName(f"scroll{title}")
        scroll.setWidgetResizable(True)
        scroll_policy = getattr(self.qt.core.Qt, "ScrollBarPolicy", self.qt.core.Qt)
        scroll.setHorizontalScrollBarPolicy(scroll_policy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(scroll_policy.ScrollBarAsNeeded)
        scroll.setFrameShape(getattr(self.w.QFrame, "Shape", self.w.QFrame).NoFrame)
        canvas = self.w.QWidget(scroll)
        canvas.setObjectName(f"{title}PageContent")
        canvas.setProperty("pageCanvas", True)
        inner = self.w.QVBoxLayout(canvas)
        inner.setContentsMargins(18, 18, 18, 24)
        inner.setSpacing(18)
        inner.addWidget(_label(self.w, PAGE_LABELS.get(title, title), "pageTitle"))
        purpose_label = _label(self.w, purpose, "pagePurpose")
        purpose_label.setWordWrap(True)
        inner.addWidget(purpose_label)
        content = self.w.QVBoxLayout()
        content.setSpacing(16)
        inner.addLayout(content)
        inner.addStretch(1)
        scroll.setWidget(canvas)
        outer.addWidget(scroll)
        self.pages.addWidget(page)
        return page, content

    # --------------------------------------------------------------- Project

    def _wire_navigation(self) -> None:
        self.nav.currentRowChanged.connect(self.pages.setCurrentIndex)


__all__ = ["ShellMixin"]
